#version 330 core
#define EXPONENT 50
#define Directional 0
#define Point 1
#define Spot 2
#define Pixel 0
#define Vertex 1

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

out vec4 FragColor;

in vec3 vertex_color;
in vec3 vertex_normal;
in vec3 vertex_position;
in vec3 vertex_position2;

uniform int cur_light_mode;
uniform int fragment_mode;
uniform float shininess;
uniform vec3 light_position;
uniform vec3 light_direction;

uniform Reflection reflection;
uniform Lighting lighting;

void main() {
	// [TODO]
	// ambient
	vec3 Ia = vec3(0.15f, 0.15f, 0.15f);
	vec3 ambient = Ia * reflection.Ka;

	// diffuse
	vec3 Id = lighting.diffuse_intensity;
	vec3 normalized_N = normalize(vertex_normal);
	vec3 normalized_LD = normalize(lighting.position);
	vec3 diffuse = Id * reflection.Kd * max(dot(normalized_N, normalized_LD), 0.0);

	// attenuation
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

	vec4 color;

	if(cur_light_mode == Directional) {
		color = vec4(ambient + diffuse + specular, 0.0f);
	} else if(cur_light_mode == Point) {
		float distance = length(lighting.position - vertex_position);
		float dist = constant + linear * distance + quadratic * distance * distance;
		float attenuation = min(1.0f / dist, 1.0f);
		color = vec4(ambient + attenuation * (diffuse + specular), 0.0f);
	} else if(cur_light_mode == Spot) {
		float spotDot = dot(normalize(vertex_position - lighting.position), normalize(lighting.direction));
		if(spotDot > cos(radians(lighting.angle))) {
			float dot = pow(spotDot, EXPONENT);
			color = vec4(ambient + dot * (diffuse + specular), 0.0f);
		} else {
			color = vec4(ambient, 0.0f);
		}
	}

	if(fragment_mode == Vertex) {
		FragColor = vec4(vertex_color, 0.0f);
	} else if(fragment_mode == Pixel) {
		FragColor = color;
	}
}
