#version 330
#define EXPONENT 50
#define Directional 0
#define Point 1
#define Spot 2

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aColor;
layout (location = 2) in vec3 aNormal;
layout (location = 3) in vec2 aTexCoord;

struct Reflection {
	vec3 Ka;
	vec3 Kd;
	vec3 Ks;
};

struct Lighting {
	vec3 position;
	vec3 direction;
	vec3 diffuse_intensity;
	float angle;
};

struct Matrix {
	mat4 m;
	mat4 v;
	mat4 p;
};

out vec2 texCoord;
out vec3 vertex_color;
out vec3 vertex_normal;
out vec3 vertex_position;
out vec3 vertex_position2;

uniform int cur_light_mode;
uniform float shininess;
uniform int isEye;
uniform float textureX;
uniform float textureY;

uniform Reflection reflection;
uniform Lighting lighting;
uniform Matrix matrix;

// [TODO] passing uniform variable for texture coordinate offset

void main() 
{
	// [TODO]
	if(isEye == 1) {
		texCoord = aTexCoord + vec2(textureX / 2, textureY / (-4));
	} else {
		texCoord = aTexCoord;
	}
	gl_Position = matrix.p * matrix.v * matrix.m * vec4(aPos, 1.0);
	vertex_position = (matrix.m * vec4(aPos.x, aPos.y, aPos.z, 1.0)).xyz;
	vertex_position2 = (matrix.v * matrix.m * vec4(aPos.x, aPos.y, aPos.z, 1.0)).xyz;
	vertex_normal = (transpose(inverse(matrix.v * matrix.m)) * vec4(aNormal, 0.0)).xyz;

	//ambient
	vec3 Ia = vec3(0.15f, 0.15f, 0.15f);
	vec3 ambient = Ia * reflection.Ka;

	//diffuse
	vec3 Id = lighting.diffuse_intensity;
	vec3 normalized_N = normalize(vertex_normal);
	vec3 normalized_LD = normalize(lighting.position);
	vec3 diffuse = Id * reflection.Kd * max(dot(normalized_N, normalized_LD), 0.0);

	//attenuation
	float constant;
	float linear;
	float quadratic;
	
	if(cur_light_mode == Point) {
		constant = 0.01;
		linear = 0.8;
		quadratic = 0.1;
	} else if(cur_light_mode == Spot) {
		constant = 0.05;
		linear = 0.3;
		quadratic = 0.6;
	}

	// specular
	vec3 Is = vec3(1.0f, 1.0f, 1.0f);
	vec3 normalized_Rp = normalize(lighting.position - vertex_position2);
	vec3 specular = Is * reflection.Ks * pow(max(dot(normalized_Rp, normalized_N), 0.0), shininess);

	if(cur_light_mode == Directional) {
		vertex_color = ambient + diffuse + specular;
	} else if(cur_light_mode == Point) {
		float distance = length(lighting.position - vertex_position);
		float dist = constant + linear * distance + quadratic * distance * distance;
		float attenuation = min(1.0f / dist, 1.0f);
		vertex_color = ambient + attenuation * (diffuse + specular);
	} else if(cur_light_mode == Spot) {
		float spotDot = dot(normalize(vertex_position - lighting.position), normalize(lighting.direction));
		if(spotDot > cos(radians(lighting.angle))) {
			vertex_color = ambient + pow(spotDot, EXPONENT) * (diffuse + specular);
		} else {
			vertex_color = ambient;
		}
	}
}
