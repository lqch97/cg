#include <iostream> 
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "textfile.h"

#include "Vectors.h"
#include "Matrices.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

#ifndef max
# define max(a,b) (((a)>(b))?(a):(b))
# define min(a,b) (((a)<(b))?(a):(b))
#endif

using namespace std;

// Default window size
const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 800;

// window size in pixel
int window_width_in_pixel;
int window_height_in_pixel;

// Basic unit for angle degree
const float ANGLE_DEGREE = 3.14159265359 / 180;

// variables for mouse control
bool mouse_pressed = false;
int starting_press_x = -1;
int starting_press_y = -1;

enum TransMode
{
	GeoTranslation = 0,
	GeoRotation = 1,
	GeoScaling = 2,
	LightEdit = 3,
	ShininessEdit = 4,
};

struct Uniform
{
    // flags to control the flow of shading
    GLint iLocLightMode;
    GLint iLocShadingMode;
    
    // matirx
	GLint iLocMVP;
    GLint iLocMV;
    GLint iLocV;
    GLint iLocNormTrnas;
    
    // lighting properties (general)
    GLint iLocPosition;
    GLint iLocDiffuseIntensity;
    GLint iLocAmbientIntensity;
    GLint iLocSpecularIntensity;
    
    // attenuation
    GLint iLocConstant;
    GLint iLocLinear;
    GLint iLocQuadratic;
    
    // lighting properties (spot light)
    GLint iLocDirection;
    GLint iLocExponent;
    GLint iLocCutoff;
    
    // material properties
    GLint iLocKa;
    GLint iLocKd;
    GLint iLocKs;
    GLint iLocShininess;
    
    // [my TODO] for debug, delete later
    GLint iLocFlag;
};
Uniform uniform;

// properties for light source
struct Light {
    // general
    Vector3 position;
    Vector3 diffuseIntensity;
    Vector3 ambientIntensity;
    Vector3 specularIntensity;
};

struct DirectLight: Light {}; // no new property
DirectLight DL;

struct PositionLight: Light {
    // attenuation
    float constant;
    float linear;
    float quadratic;
};
PositionLight PL;

struct SpotLight: PositionLight {
    // for spotlight
    Vector3 direction;
    float exponent;
    float cutoff;
};
SpotLight SL;

vector<string> filenames; // .obj filename list

struct PhongMaterial
{
	Vector3 Ka;
	Vector3 Kd;
	Vector3 Ks;
};

typedef struct
{
	GLuint vao;
	GLuint vbo;
	GLuint vboTex;
	GLuint ebo;
	GLuint p_color;
	int vertex_count;
	GLuint p_normal;
	PhongMaterial material;
	int indexCount;
	GLuint m_texture;
} Shape;

struct model
{
	Vector3 position = Vector3(0, 0, 0);
	Vector3 scale = Vector3(1, 1, 1);
	Vector3 rotation = Vector3(0, 0, 0);	// Euler form

	vector<Shape> shapes;
};
vector<model> models;

struct camera
{
	Vector3 position;
	Vector3 center;
	Vector3 up_vector;
};
camera main_camera;

struct project_setting
{
	GLfloat nearClip, farClip;
	GLfloat fovy;
	GLfloat aspect;
	GLfloat left, right, top, bottom;
};
project_setting proj;

TransMode cur_trans_mode = GeoTranslation;
int cur_light_mode = 1; // [0, 1, 2] = [direct, point, spot], the light mode be used now [myTODO] need a better representation
int cur_idx = 0; // represent which model should be rendered now

Matrix4 view_matrix;
Matrix4 project_matrix;

// material property for specular
GLfloat shininess;

// [my TODO] for debuggin, delete later
int flag = 3;


// [DO] given a translation vector then output a Matrix4 (Translation Matrix)
Matrix4 translate(Vector3 vec)
{
    return Matrix4(
        1, 0, 0, vec.x,
        0, 1, 0, vec.y,
        0, 0, 1, vec.z,
        0, 0, 0, 1
    );
}

// [DO] given a scaling vector then output a Matrix4 (Scaling Matrix)
Matrix4 scaling(Vector3 vec)
{
    return Matrix4(
        vec.x, 0, 0, 0,
        0, vec.y, 0, 0,
        0, 0, vec.z, 0,
        0, 0, 0, 1
    );
}


// [DO] given a float value then ouput a rotation matrix alone axis-X (rotate alone axis-X)
Matrix4 rotateX(GLfloat val)
{
    return  Matrix4(
        1, 0, 0, 0,
        0, cos(val), -sin(val), 0,
        0, sin(val), cos(val), 0,
        0, 0, 0, 1
    );
}

// [DO] given a float value then ouput a rotation matrix alone axis-Y (rotate alone axis-Y)
Matrix4 rotateY(GLfloat val)
{
    return Matrix4(
        cos(val), 0, sin(val), 0,
        0, 1, 0, 0,
        -sin(val), 0, cos(val), 0,
        0, 0, 0, 1
    );
}

// [DO] given a float value then ouput a rotation matrix alone axis-Z (rotate alone axis-Z)
Matrix4 rotateZ(GLfloat val)
{
    return Matrix4(
        cos(val), -sin(val), 0, 0,
        sin(val), cos(val), 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    );
}

Matrix4 rotate(Vector3 vec)
{
	return rotateX(vec.x)*rotateY(vec.y)*rotateZ(vec.z);
}

// [DO] compute viewing matrix accroding to the setting of main_camera
void setViewingMatrix()
{
    Vector3 p1_p2 = main_camera.center - main_camera.position,
                p1_p3 = main_camera.up_vector,
                new_z = -(p1_p2).normalize(),
                new_x = p1_p2.cross(p1_p3).normalize(),
                new_y = new_z.cross(new_x);

    Matrix4 T = translate(-main_camera.position);
    Matrix4 R = Matrix4(
        new_x[0], new_x[1], new_x[2], 0,
        new_y[0], new_y[1], new_y[2], 0,
        new_z[0], new_z[1], new_z[2], 0,
        0, 0, 0, 1
    );

    view_matrix = R * T;
}

// [DO] compute persepective projection matrix
void setPerspective()
{
    float f = 1 / (tan(proj.fovy / 2 * ANGLE_DEGREE)); // f = cot(fovy / 2)
    project_matrix = Matrix4(
        f / proj.aspect, 0, 0, 0,
        0, f, 0, 0,
        0, 0, (proj.farClip + proj.nearClip) / (proj.nearClip - proj.farClip), (2 * proj.farClip * proj.nearClip) / (proj.nearClip - proj.farClip),
        0, 0, -1, 0
    );
}

void setGLMatrix(GLfloat* glm, Matrix4& m) {
	glm[0] = m[0];  glm[4] = m[1];   glm[8] = m[2];    glm[12] = m[3];
	glm[1] = m[4];  glm[5] = m[5];   glm[9] = m[6];    glm[13] = m[7];
	glm[2] = m[8];  glm[6] = m[9];   glm[10] = m[10];   glm[14] = m[11];
	glm[3] = m[12];  glm[7] = m[13];  glm[11] = m[14];   glm[15] = m[15];
}

// Call back function for window reshape
void ChangeSize(GLFWwindow* window, int width, int height)
{
	// [DO] change your aspect ratio
    proj.aspect = ((float)window_width_in_pixel / 2) / (float)window_height_in_pixel: // [my TODO] width of aspect ratio
    setPerspective();

    // reset window size
    window_width_in_pixel = width;
    window_height_in_pixel = height;
}

// set properties to uniform variable in shader
void setUniforms() {
    // [DO] update translation, rotation and scaling
    // [DO] multiply all the matrix
    model cur_model = models[cur_idx];
    Matrix4 T = translate(cur_model.position),
            R = rotate(cur_model.rotation),
            S = scaling(cur_model.scale);
    
    // temp variable for passing matrix to shader
    GLfloat temp[16];
    
    // set MV matrix
    Matrix4 MV =  view_matrix * T * R * S;
    setGLMatrix(temp, MV); // row-major ---> column-major
    glUniformMatrix4fv(uniform.iLocMV, 1, GL_FALSE, temp); // use uniform to send mvp to vertex shader
    
    // set viewing matrix
    setGLMatrix(temp, view_matrix); // row-major ---> column-major
    glUniformMatrix4fv(uniform.iLocV, 1, GL_FALSE, temp);
    
    // set MVP matrix
    Matrix4 MVP = project_matrix * MV;
    setGLMatrix(temp, MVP); // row-major ---> column-major
    glUniformMatrix4fv(uniform.iLocMVP, 1, GL_FALSE, temp);
    
    // set mormal transformation matrix
    Matrix4 NORM_TRANS = MV.invert().transpose();
    setGLMatrix(temp, NORM_TRANS); // row-major ---> column-major
    glUniformMatrix4fv(uniform.iLocNormTrnas, 1, GL_FALSE, temp);
    
    switch (cur_light_mode) {
        case 0: //directLight
            // general
            glUniform3f(uniform.iLocPosition, DL.position.x, DL.position.y, DL.position.z);
            glUniform3f(uniform.iLocAmbientIntensity, DL.ambientIntensity.x, DL.ambientIntensity.y, DL.ambientIntensity.z);
            glUniform3f(uniform.iLocDiffuseIntensity, DL.diffuseIntensity.x, DL.diffuseIntensity.y, DL.diffuseIntensity.z);
            glUniform3f(uniform.iLocSpecularIntensity, DL.specularIntensity.x, DL.specularIntensity.y, DL.specularIntensity.z);
            break;
        
        case 1: //positionLight
            // general
            glUniform3f(uniform.iLocPosition, PL.position.x, PL.position.y, PL.position.z);
            glUniform3f(uniform.iLocAmbientIntensity, PL.ambientIntensity.x, PL.ambientIntensity.y, PL.ambientIntensity.z);
            glUniform3f(uniform.iLocDiffuseIntensity, PL.diffuseIntensity.x, PL.diffuseIntensity.y, PL.diffuseIntensity.z);
            glUniform3f(uniform.iLocSpecularIntensity, PL.specularIntensity.x, PL.specularIntensity.y, PL.specularIntensity.z);
            // attenuation
            glUniform1f(uniform.iLocConstant, PL.constant);
            glUniform1f(uniform.iLocLinear, PL.linear);
            glUniform1f(uniform.iLocQuadratic, PL.quadratic);
            break;
            
        case 2: // spotLight
            // general
            glUniform3f(uniform.iLocPosition, SL.position.x, SL.position.y, SL.position.z);
            glUniform3f(uniform.iLocAmbientIntensity, SL.ambientIntensity.x, SL.ambientIntensity.y, SL.ambientIntensity.z);
            glUniform3f(uniform.iLocDiffuseIntensity, SL.diffuseIntensity.x, SL.diffuseIntensity.y, SL.diffuseIntensity.z);
            glUniform3f(uniform.iLocSpecularIntensity, SL.specularIntensity.x, SL.specularIntensity.y, SL.specularIntensity.z);
            // attenuation
            glUniform1f(uniform.iLocConstant, SL.constant);
            glUniform1f(uniform.iLocLinear, SL.linear);
            glUniform1f(uniform.iLocQuadratic, SL.quadratic);
            // spotlight specific setting
            glUniform3f(uniform.iLocDirection, SL.direction.x, SL.direction.y, SL.direction.z);
            glUniform1f(uniform.iLocExponent, SL.exponent);
            glUniform1f(uniform.iLocCutoff, SL.cutoff); 
            break;
    }
    
    // set light mode
    glUniform1i(uniform.iLocLightMode, cur_light_mode);
    
    // [my TODO] debug
    glUniform1i(uniform.iLocFlag, flag);
}

// Render function for display rendering
void RenderScene(void) {	
	// clear canvas
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    
    setUniforms();
    
    for (int i = 0; i < models[cur_idx].shapes.size(); i++) {
        // material properties
        glUniform3fv(uniform.iLocKa, 1, reinterpret_cast<GLfloat*>(&models[cur_idx].shapes[i].material.Ka[0]));
        glUniform3fv(uniform.iLocKd, 1, reinterpret_cast<GLfloat*>(&models[cur_idx].shapes[i].material.Kd[0]));
        glUniform3fv(uniform.iLocKs, 1, reinterpret_cast<GLfloat*>(&models[cur_idx].shapes[i].material.Ks[0]));
        glUniform1f(uniform.iLocShininess, shininess);
                
        // draw left hand side viewport in vertex lighting
        glUniform1i(uniform.iLocShadingMode, 0); // set shading mode, 0 => vertex shading
        glViewport(0, 0, (GLsizei)(window_width_in_pixel / 2), window_height_in_pixel);
        glBindVertexArray(models[cur_idx].shapes[i].vao);
        glDrawArrays(GL_TRIANGLES, 0, models[cur_idx].shapes[i].vertex_count);
                
        // draw right hand side viewport in pixel lighting
        glUniform1i(uniform.iLocShadingMode, 1); // set shading mode, 1 => fragment shading
        glViewport((GLsizei)(window_width_in_pixel / 2), 0, (GLsizei)(window_width_in_pixel / 2), window_height_in_pixel);
        glBindVertexArray(models[cur_idx].shapes[i].vao);
        glDrawArrays(GL_TRIANGLES, 0, models[cur_idx].shapes[i].vertex_count);
    }
}


void KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	// [DO] Call back function for keyboard
    if( !(action == GLFW_PRESS) ) return;
        
    unsigned int num_models = (unsigned int) models.size();
        
    switch (key) {
        case GLFW_KEY_ESCAPE:
            exit(0);
                    
        case GLFW_KEY_Z:
            cur_idx = (cur_idx + num_models - 1) % num_models;
            break;
                    
        case GLFW_KEY_X:
            cur_idx = (cur_idx + 1) % num_models;
            break;
            
        case GLFW_KEY_T:
            cur_trans_mode = GeoTranslation;
            break;
                                                                                
        case GLFW_KEY_S:
            cur_trans_mode = GeoScaling;
            break;
            
        case GLFW_KEY_R:
            cur_trans_mode = GeoRotation;
            break;
                                                                                
        case GLFW_KEY_L:
            cur_light_mode = (cur_light_mode + 1) % 3;
            cout << "Light Mode: "
                 << ((cur_light_mode == 0) ? "Directional" :
                     (cur_light_mode == 1) ? "Positional" :
                                             "Spot")
                 << " Light" << endl;
            break;
                                                                                
        case GLFW_KEY_K:
            cur_trans_mode = LightEdit;
            break;
                                                                                
        case GLFW_KEY_J:
            cur_trans_mode = ShininessEdit;
            break;
            
        // [my TODO] debug   
        case GLFW_KEY_F:
            flag = (flag + 1) % 4;
            break;
            
    }
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	// [DO] scroll up positive, otherwise it would be negtive
    const float translation_factor = 0.01;
    const float scaling_factor = 0.01;
    const float rotation_factor = 0.2;
    const float shininess_changing_factor = 2.0;
    const float cutoff_changing_factor = 0.5;
    const float diffuse_changing_factor = 0.1;
    
    switch (cur_trans_mode) {
        case GeoTranslation:
            models[cur_idx].position.z += yoffset * translation_factor;
            break;
                        
        case GeoScaling:
            models[cur_idx].scale.z += yoffset * scaling_factor;
            break;
                        
        case GeoRotation:
            models[cur_idx].rotation.z += ANGLE_DEGREE * yoffset * rotation_factor;
            break;
            
        case ShininessEdit:
            shininess = max(shininess + yoffset * shininess_changing_factor, 1);
            break;
            
        case LightEdit:
            if(cur_light_mode == 2) { // spotlight mode
                SL.cutoff +=  yoffset * cutoff_changing_factor * ANGLE_DEGREE; // [my TODO] min, max degree ?
            } else {
                auto& diffuse = (cur_light_mode == 0) ? DL.diffuseIntensity : PL.diffuseIntensity;
                diffuse += Vector3(yoffset * diffuse_changing_factor, yoffset * diffuse_changing_factor, yoffset * diffuse_changing_factor);
            }
            break;
    }
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	// [DO] mouse press callback function
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            mouse_pressed = true;
        }
        if (action == GLFW_RELEASE) {
            mouse_pressed = false;
            starting_press_x = -1;
            starting_press_y = -1;
        }
    }
}

static void cursor_pos_callback(GLFWwindow* window, double xpos, double ypos)
{
	// [DO] cursor position callback function
    if(!mouse_pressed) return;
        
    if(starting_press_x == -1 || starting_press_y == -1) {
        starting_press_x = xpos;
        starting_press_y = ypos;
        return;
    }
        
    int delta_x = xpos - starting_press_x;
    int delta_y = -(ypos - starting_press_y); // move up: - -> +
        
    const float translation_factor = 0.01;
    const float scaling_factor = 0.01;
    const float rotation_factor = 0.2;
    const float light_translation_factor = 0.01;
    
    switch (cur_trans_mode) {
        case GeoTranslation:
            models[cur_idx].position.x += delta_x * translation_factor;
            models[cur_idx].position.y += delta_y * translation_factor;
            break;
                        
        case GeoScaling:
            models[cur_idx].scale.x -= delta_x * scaling_factor;
            models[cur_idx].scale.y += delta_y * scaling_factor;
            break;
                        
        case GeoRotation:
            models[cur_idx].rotation.x += ANGLE_DEGREE * delta_y * rotation_factor;
            models[cur_idx].rotation.y -= ANGLE_DEGREE * delta_x * rotation_factor;
            break;
            
        case LightEdit:
            auto& position = (cur_light_mode == 0) ? DL.position:
                             (cur_light_mode == 1) ? PL.position:
                                                     SL.position;
            position.x += delta_x * light_translation_factor;
            position.y += delta_y * light_translation_factor;
            break;
    }
    
    starting_press_x = xpos;
    starting_press_y = ypos;
}

void setILoc(GLint program) {
    
    // flags to control the flow of shading
    uniform.iLocLightMode = glGetUniformLocation(program, "lightMode");
    uniform.iLocShadingMode = glGetUniformLocation(program, "shadingMode");
    
    // matrix
    uniform.iLocMVP = glGetUniformLocation(program, "mvp");
    uniform.iLocMV = glGetUniformLocation(program, "mv");
    uniform.iLocV = glGetUniformLocation(program, "v");
    uniform.iLocNormTrnas = glGetUniformLocation(program, "normTrans");

    // general lighting properties
    uniform.iLocPosition = glGetUniformLocation(program, "position");
    uniform.iLocDiffuseIntensity = glGetUniformLocation(program, "diffuseIntensity");
    uniform.iLocAmbientIntensity = glGetUniformLocation(program, "ambientIntensity");
    uniform.iLocSpecularIntensity = glGetUniformLocation(program, "specularIntensity");

    // attenuation
    uniform.iLocConstant = glGetUniformLocation(program, "constant"); // shiniess ?
    uniform.iLocLinear = glGetUniformLocation(program, "linear");
    uniform.iLocQuadratic = glGetUniformLocation(program, "quadratic");

    // for spotlight
    uniform.iLocDirection = glGetUniformLocation(program, "direction");
    uniform.iLocExponent = glGetUniformLocation(program, "exponent");
    uniform.iLocCutoff = glGetUniformLocation(program, "cutoff");
    
    // material properties
    uniform.iLocKa = glGetUniformLocation(program, "Ka");
    uniform.iLocKd = glGetUniformLocation(program, "Kd");
    uniform.iLocKs = glGetUniformLocation(program, "Ks");
    uniform.iLocShininess = glGetUniformLocation(program, "shininess");
    
    // var for debug [my TODO]
    uniform.iLocFlag = glGetUniformLocation(program, "flag");
}

void setShaders()
{
	GLuint v, f, p;
	char *vs = NULL;
	char *fs = NULL;

	v = glCreateShader(GL_VERTEX_SHADER);
	f = glCreateShader(GL_FRAGMENT_SHADER);

	vs = textFileRead("shader.vs");
	fs = textFileRead("shader.fs");

	glShaderSource(v, 1, (const GLchar**)&vs, NULL);
	glShaderSource(f, 1, (const GLchar**)&fs, NULL);

	free(vs);
	free(fs);

	GLint success;
	char infoLog[1000];
	// compile vertex shader
	glCompileShader(v);
	// check for shader compile errors
	glGetShaderiv(v, GL_COMPILE_STATUS, &success);
	if (!success)
	{
		glGetShaderInfoLog(v, 1000, NULL, infoLog);
		std::cout << "ERROR: VERTEX SHADER COMPILATION FAILED\n" << infoLog << std::endl;
	}

	// compile fragment shader
	glCompileShader(f);
	// check for shader compile errors
	glGetShaderiv(f, GL_COMPILE_STATUS, &success);
	if (!success)
	{
		glGetShaderInfoLog(f, 1000, NULL, infoLog);
		std::cout << "ERROR: FRAGMENT SHADER COMPILATION FAILED\n" << infoLog << std::endl;
	}

	// create program object
	p = glCreateProgram();

	// attach shaders to program object
	glAttachShader(p,f);
	glAttachShader(p,v);

	// link program
	glLinkProgram(p);
	// check for linking errors
	glGetProgramiv(p, GL_LINK_STATUS, &success);
	if (!success) {
		glGetProgramInfoLog(p, 1000, NULL, infoLog);
		std::cout << "ERROR: SHADER PROGRAM LINKING FAILED\n" << infoLog << std::endl;
	}

	glDeleteShader(v);
	glDeleteShader(f);

    setILoc(p);

	if (success)
		glUseProgram(p);
    else
    {
        system("pause");
        exit(123);
    }
}

void normalization(tinyobj::attrib_t* attrib, vector<GLfloat>& vertices, vector<GLfloat>& colors, vector<GLfloat>& normals, tinyobj::shape_t* shape)
{
	vector<float> xVector, yVector, zVector;
	float minX = 10000, maxX = -10000, minY = 10000, maxY = -10000, minZ = 10000, maxZ = -10000;

	// find out min and max value of X, Y and Z axis
	for (int i = 0; i < attrib->vertices.size(); i++)
	{
		//maxs = max(maxs, attrib->vertices.at(i));
		if (i % 3 == 0)
		{

			xVector.push_back(attrib->vertices.at(i));

			if (attrib->vertices.at(i) < minX)
			{
				minX = attrib->vertices.at(i);
			}

			if (attrib->vertices.at(i) > maxX)
			{
				maxX = attrib->vertices.at(i);
			}
		}
		else if (i % 3 == 1)
		{
			yVector.push_back(attrib->vertices.at(i));

			if (attrib->vertices.at(i) < minY)
			{
				minY = attrib->vertices.at(i);
			}

			if (attrib->vertices.at(i) > maxY)
			{
				maxY = attrib->vertices.at(i);
			}
		}
		else if (i % 3 == 2)
		{
			zVector.push_back(attrib->vertices.at(i));

			if (attrib->vertices.at(i) < minZ)
			{
				minZ = attrib->vertices.at(i);
			}

			if (attrib->vertices.at(i) > maxZ)
			{
				maxZ = attrib->vertices.at(i);
			}
		}
	}

	float offsetX = (maxX + minX) / 2;
	float offsetY = (maxY + minY) / 2;
	float offsetZ = (maxZ + minZ) / 2;

	for (int i = 0; i < attrib->vertices.size(); i++)
	{
		if (offsetX != 0 && i % 3 == 0)
		{
			attrib->vertices.at(i) = attrib->vertices.at(i) - offsetX;
		}
		else if (offsetY != 0 && i % 3 == 1)
		{
			attrib->vertices.at(i) = attrib->vertices.at(i) - offsetY;
		}
		else if (offsetZ != 0 && i % 3 == 2)
		{
			attrib->vertices.at(i) = attrib->vertices.at(i) - offsetZ;
		}
	}

	float greatestAxis = maxX - minX;
	float distanceOfYAxis = maxY - minY;
	float distanceOfZAxis = maxZ - minZ;

	if (distanceOfYAxis > greatestAxis)
	{
		greatestAxis = distanceOfYAxis;
	}

	if (distanceOfZAxis > greatestAxis)
	{
		greatestAxis = distanceOfZAxis;
	}

	float scale = greatestAxis / 2;

	for (int i = 0; i < attrib->vertices.size(); i++)
	{
		//std::cout << i << " = " << (double)(attrib.vertices.at(i) / greatestAxis) << std::endl;
		attrib->vertices.at(i) = attrib->vertices.at(i) / scale;
	}
	size_t index_offset = 0;
	for (size_t f = 0; f < shape->mesh.num_face_vertices.size(); f++) {
		int fv = shape->mesh.num_face_vertices[f];

		// Loop over vertices in the face.
		for (size_t v = 0; v < fv; v++) {
			// access to vertex
			tinyobj::index_t idx = shape->mesh.indices[index_offset + v];
			vertices.push_back(attrib->vertices[3 * idx.vertex_index + 0]);
			vertices.push_back(attrib->vertices[3 * idx.vertex_index + 1]);
			vertices.push_back(attrib->vertices[3 * idx.vertex_index + 2]);
			// Optional: vertex colors
			colors.push_back(attrib->colors[3 * idx.vertex_index + 0]);
			colors.push_back(attrib->colors[3 * idx.vertex_index + 1]);
			colors.push_back(attrib->colors[3 * idx.vertex_index + 2]);
			// Optional: vertex normals
			if (idx.normal_index >= 0) {
				normals.push_back(attrib->normals[3 * idx.normal_index + 0]);
				normals.push_back(attrib->normals[3 * idx.normal_index + 1]);
				normals.push_back(attrib->normals[3 * idx.normal_index + 2]);
			}
		}
		index_offset += fv;
	}
}

string GetBaseDir(const string& filepath) {
	if (filepath.find_last_of("/\\") != std::string::npos)
		return filepath.substr(0, filepath.find_last_of("/\\"));
	return "";
}

void LoadModels(string model_path)
{
	vector<tinyobj::shape_t> shapes;
	vector<tinyobj::material_t> materials;
	tinyobj::attrib_t attrib;
	vector<GLfloat> vertices;
	vector<GLfloat> colors;
	vector<GLfloat> normals;

	string err;
	string warn;

	string base_dir = GetBaseDir(model_path); // handle .mtl with relative path

#ifdef _WIN32
	base_dir += "\\";
#else
	base_dir += "/";
#endif

	bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, model_path.c_str(), base_dir.c_str());

	if (!warn.empty()) {
		cout << warn << std::endl;
	}

	if (!err.empty()) {
		cerr << err << std::endl;
	}

	if (!ret) {
		exit(1);
	}

	printf("Load Models Success ! Shapes size %d Material size %d\n", int(shapes.size()), int(materials.size()));
	model tmp_model;

	vector<PhongMaterial> allMaterial;
	for (int i = 0; i < materials.size(); i++)
	{
		PhongMaterial material;
		material.Ka = Vector3(materials[i].ambient[0], materials[i].ambient[1], materials[i].ambient[2]);
		material.Kd = Vector3(materials[i].diffuse[0], materials[i].diffuse[1], materials[i].diffuse[2]);
		material.Ks = Vector3(materials[i].specular[0], materials[i].specular[1], materials[i].specular[2]);
		allMaterial.push_back(material);
	}

	for (int i = 0; i < shapes.size(); i++)
	{

		vertices.clear();
		colors.clear();
		normals.clear();
		normalization(&attrib, vertices, colors, normals, &shapes[i]);
		// printf("Vertices size: %d", vertices.size() / 3);

		Shape tmp_shape;
		glGenVertexArrays(1, &tmp_shape.vao);
		glBindVertexArray(tmp_shape.vao);

		glGenBuffers(1, &tmp_shape.vbo);
		glBindBuffer(GL_ARRAY_BUFFER, tmp_shape.vbo);
		glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(GL_FLOAT), &vertices.at(0), GL_STATIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		tmp_shape.vertex_count = vertices.size() / 3;

		glGenBuffers(1, &tmp_shape.p_color);
		glBindBuffer(GL_ARRAY_BUFFER, tmp_shape.p_color);
		glBufferData(GL_ARRAY_BUFFER, colors.size() * sizeof(GL_FLOAT), &colors.at(0), GL_STATIC_DRAW);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

		glGenBuffers(1, &tmp_shape.p_normal);
		glBindBuffer(GL_ARRAY_BUFFER, tmp_shape.p_normal);
		glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(GL_FLOAT), &normals.at(0), GL_STATIC_DRAW);
		glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, 0);

		glEnableVertexAttribArray(0);
		glEnableVertexAttribArray(1);
		glEnableVertexAttribArray(2);

		// not support per face material, use material of first face
		if (allMaterial.size() > 0)
			tmp_shape.material = allMaterial[shapes[i].mesh.material_ids[0]];
		tmp_model.shapes.push_back(tmp_shape);
	}
	shapes.clear();
	materials.clear();
	models.push_back(tmp_model);
}

void initParameter()
{
	// [DO] Setup some parameters if you need
	proj.left = -1;
	proj.right = 1;
	proj.top = 1;
	proj.bottom = -1;
	proj.nearClip = 0.001;
	proj.farClip = 100.0;
	proj.fovy = 80;
	proj.aspect = (float)window_width_in_pixel / (float)window_height_in_pixel;

	main_camera.position = Vector3(0.0f, 0.0f, 2.0f);
	main_camera.center = Vector3(0.0f, 0.0f, 0.0f);
	main_camera.up_vector = Vector3(0.0f, 1.0f, 0.0f);

	setViewingMatrix();
	setPerspective();	//set default projection matrix as perspective matrix
    
    // properties for lighting
    shininess = 64;

    DL.position = Vector3(1.0, 1.0, 1.0);
    DL.ambientIntensity = Vector3(0.15, 0.15, 0.15);
    DL.diffuseIntensity = Vector3(1.0, 1.0, 1.0);
    DL.specularIntensity = Vector3(1.0, 1.0, 1.0);

    PL.position = Vector3(0.0, 2.0, 1.0);
    PL.ambientIntensity = Vector3(0.15, 0.15, 0.15);
    PL.diffuseIntensity = Vector3(1.0, 1.0, 1.0);
    PL.specularIntensity = Vector3(1.0, 1.0, 1.0);
    PL.constant = 0.01;
    PL.linear = 0.8;
    PL.quadratic = 0.1;

    SL.position = Vector3(0.0, 0.0, 2.0);
    SL.ambientIntensity = Vector3(0.15, 0.15, 0.15);
    SL.diffuseIntensity = Vector3(1.0, 1.0, 1.0);
    SL.specularIntensity = Vector3(1.0, 1.0, 1.0);
    SL.constant = 0.05;
    SL.linear = 0.3;
    SL.quadratic = 0.6;
    SL.direction = Vector3(0.0, 0.0, -1.0);
    SL.exponent = 50.0;
    SL.cutoff = 30 * ANGLE_DEGREE;
}

void setupRC()
{
    // setup shaders
    setShaders();
    initParameter();

    // OpenGL States and Values
    glClearColor(0.2, 0.2, 0.2, 1.0);
    vector<string> model_list{ "../NormalModels/bunny5KN.obj", "../NormalModels/dragon10KN.obj", "../NormalModels/lucy25KN.obj", "../NormalModels/teapot4KN.obj", "../NormalModels/dolphinN.obj"};
    // [DO] Load five model at here
    for(auto m: model_list) {
        LoadModels(m);
    }
}

void glPrintContextInfo(bool printExtension)
{
	cout << "GL_VENDOR = " << (const char*)glGetString(GL_VENDOR) << endl;
	cout << "GL_RENDERER = " << (const char*)glGetString(GL_RENDERER) << endl;
	cout << "GL_VERSION = " << (const char*)glGetString(GL_VERSION) << endl;
	cout << "GL_SHADING_LANGUAGE_VERSION = " << (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;
	if (printExtension)
	{
		GLint numExt;
		glGetIntegerv(GL_NUM_EXTENSIONS, &numExt);
		cout << "GL_EXTENSIONS =" << endl;
		for (GLint i = 0; i < numExt; i++)
		{
			cout << "\t" << (const char*)glGetStringi(GL_EXTENSIONS, i) << endl;
		}
	}
}


int main(int argc, char **argv)
{
    // initial glfw
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // fix compilation on OS X
#endif

    
    // create window
	GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "110062653 HW2", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    
    // set window size in pixel
    glfwGetFramebufferSize(window, &window_width_in_pixel, &window_height_in_pixel);

    
    // load OpenGL function pointer
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }
    
	// register glfw callback functions
    glfwSetKeyCallback(window, KeyCallback);
	glfwSetScrollCallback(window, scroll_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	glfwSetCursorPosCallback(window, cursor_pos_callback);

    glfwSetFramebufferSizeCallback(window, ChangeSize);
	glEnable(GL_DEPTH_TEST);
	// Setup render context
	setupRC();

	// main loop
    while (!glfwWindowShouldClose(window))
    {
        // render
        RenderScene();
        
        // swap buffer from back to front
        glfwSwapBuffers(window);
        
        // Poll input event
        glfwPollEvents();
    }
	
	// just for compatibiliy purposes
	return 0;
}
