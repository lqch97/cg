#version 330 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aColor;
layout (location = 2) in vec3 aNormal;

// output datas
out vec3 vertex_color;
out vec3 vertex_normal;

// transformation matrix
uniform mat4 mvp;
uniform mat4 mv;
uniform mat4 v; // viewing matrix
uniform mat4 normTrans;

// light properites (global)
uniform vec3 viewPos; // [my TODO] need to be deleted
uniform int lightMode;

// lighting properties (general)
uniform vec3 position;
uniform vec3 diffuseIntensity;
uniform vec3 ambientIntensity;
uniform vec3 specularIntensity;

// attenuation
uniform float constant;
uniform float linear;
uniform float quadratic;

// lighting properties (spot light)
uniform vec3 direction;
uniform float exponent; // [my TODO] int for float ?
uniform float cutoff;

// material properties
uniform vec3 Ka;
uniform vec3 Kd;
uniform vec3 Ks;
uniform vec3 shininess;


void main()
{
	// [TODO]
    // MVP transformation on positions
	gl_Position = mvp * vec4(aPos, 1.0);
    
    // normal transformation`
    vertex_normal = normalize( (normTrans * vec4(aNormal, 1.0)).xyz );
    
    // calculate light_position, viewing_position, vertex_position
    vec3 light_pos = (v * vec4(position, 1.0)).xyz;
    vec3 vertex_pos = (mv * vec4(aPos, 1.0)).xyz;
    vec3 view_pos = vec3(0, 0, 0); // because we are in viewing space
    
    // calculate diffuse
    vec3 light_vector = normalize( light_pos - vertex_pos );
    float diffuse = max( dot(light_vector, vertex_normal), 0 );
    
    // calculate specular
    vec3 view_vector = normalize( view_pos - vertex_pos );
    vec3 halfway_vector = normalize( light_vector + view_vector );
    float specular = max( dot(halfway_vector, vertex_normal), 0 );
    
    // attenuation
    float dis = length(light_vector); // distance
    float attenuation = (lightMode == 0) ? // directional light or not
                        1 / (constant+ linear * dis + quadratic * dis * dis) : 1;
    
    
    vec3 result = ambientIntensity * Ka
            + diffuse * diffuseIntensity * Kd * attenuation
            + specular * specularIntensity * Ks * attenuation;
    
//    vertex_color = result;
    vertex_color = position;
}

