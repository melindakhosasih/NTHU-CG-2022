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

int CURRENT_WIDTH = 800;
int CURRENT_HEIGHT = 800;

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

enum LightMode
{
	Directional = 0,
	Point = 1,
	Spot = 2,
};

enum FragmentMode
{
	Vertex = 1,
	Pixel = 0
};

struct Uniform
{
	GLint iLocMVP;
	GLint iLocM;
	GLint iLocV;
	GLint iLocP;
	GLint iLocKa;
	GLint iLocKd;
	GLint iLocKs;
	GLint iLocPosition;
	GLint iLocDirection;
	GLint iLocIntensity;
	GLint iLocAngle;
	GLint iLocShininess;
	GLint iLocMode;
	GLint iLocFVMode;
};
Uniform uniform;

typedef struct {
	Vector3 position;
	Vector3 direction;
	Vector3 intensity;
	GLfloat angle;
} Light;

Light lightInfo[3];
Vector4 light_position, light_direction, light_intensity;

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
GLint cur_light_mode = 0;
GLint fragment_mode = Pixel;
GLfloat shininess = 64.0f;

Matrix4 view_matrix;
Matrix4 project_matrix;

int cur_idx = 0; // represent which model should be rendered now


static GLvoid Normalize(GLfloat v[3])
{
	GLfloat l;

	l = (GLfloat)sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	v[0] /= l;
	v[1] /= l;
	v[2] /= l;
}

static GLvoid Cross(GLfloat u[3], GLfloat v[3], GLfloat n[3])
{

	n[0] = u[1] * v[2] - u[2] * v[1];
	n[1] = u[2] * v[0] - u[0] * v[2];
	n[2] = u[0] * v[1] - u[1] * v[0];
}


// [TODO] given a translation vector then output a Matrix4 (Translation Matrix)
Matrix4 translate(Vector3 vec)
{
	Matrix4 mat;

	mat = Matrix4(
		1.0f, 0.0f, 0.0f, vec.x,
		0.0f, 1.0f, 0.0f, vec.y,
		0.0f, 0.0f, 1.0f, vec.z,
		0.0f, 0.0f, 0.0f, 1.0f
	);

	return mat;
}

// [TODO] given a scaling vector then output a Matrix4 (Scaling Matrix)
Matrix4 scaling(Vector3 vec)
{
	Matrix4 mat;

	mat = Matrix4(
		vec.x, 0.0f, 0.0f, 0.0f,
		0.0f, vec.y, 0.0f, 0.0f,
		0.0f, 0.0f, vec.z, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	);

	return mat;
}


// [TODO] given a float value then ouput a rotation matrix alone axis-X (rotate alone axis-X)
Matrix4 rotateX(GLfloat val)
{
	Matrix4 mat;

	mat = Matrix4(
		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, cos(val), -sin(val), 0.0f,
		0.0f, sin(val), cos(val), 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	);

	return mat;
}

// [TODO] given a float value then ouput a rotation matrix alone axis-Y (rotate alone axis-Y)
Matrix4 rotateY(GLfloat val)
{
	Matrix4 mat;

	mat = Matrix4(
		cos(val), 0.0f, sin(val), 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		-sin(val), 0.0f, cos(val), 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	);

	return mat;
}

// [TODO] given a float value then ouput a rotation matrix alone axis-Z (rotate alone axis-Z)
Matrix4 rotateZ(GLfloat val)
{
	Matrix4 mat;

	mat = Matrix4(
		cos(val), -sin(val), 0.0f, 0.0f,
		sin(val), cos(val), 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	);

	return mat;
}

Matrix4 rotate(Vector3 vec)
{
	return rotateX(vec.x)*rotateY(vec.y)*rotateZ(vec.z);
}

// [TODO] compute viewing matrix accroding to the setting of main_camera
void setViewingMatrix()
{
	Matrix4 T, R;
	Vector3 Rx, Ry, Rz;

	T = Matrix4(
		1.0f, 0.0f, 0.0f, -main_camera.position.x,
		0.0f, 1.0f, 0.0f, -main_camera.position.y,
		0.0f, 0.0f, 1.0f, -main_camera.position.z,
		0.0f, 0.0f, 0.0f, 1.0f
	);

	Rz = main_camera.position - main_camera.center;	// camera direction
	Rz.normalize();

	Rx = main_camera.up_vector.cross(Rz);	// camera right
	Rx.normalize();

	Ry = Rz.cross(Rx);	// camera up

	R = Matrix4(
		Rx.x, Rx.y, Rx.z, 0.0f,
		Ry.x, Ry.y, Ry.z, 0.0f,
		Rz.x, Rz.y, Rz.z, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	);

	view_matrix = R * T;
}

// [TODO] compute persepective projection matrix
void setPerspective()
{
	GLfloat A, B;
	A = (proj.farClip + proj.nearClip) / (proj.nearClip - proj.farClip);
	B = 2 * proj.farClip * proj.nearClip / (proj.nearClip - proj.farClip);
	project_matrix = Matrix4(
		1 / tan(proj.fovy / 120) / proj.aspect, 0.0f, 0.0f, 0.0f,
		0.0f, 1 / tan(proj.fovy / 120), 0.0f, 0.0f,
		0.0f, 0.0f, A, B,
		0.0f, 0.0f, -1.0f, 0.0f
	);
}

void setGLMatrix(GLfloat* glm, Matrix4& m) {
	glm[0] = m[0];  glm[4] = m[1];   glm[8] = m[2];    glm[12] = m[3];
	glm[1] = m[4];  glm[5] = m[5];   glm[9] = m[6];    glm[13] = m[7];
	glm[2] = m[8];  glm[6] = m[9];   glm[10] = m[10];   glm[14] = m[11];
	glm[3] = m[12];  glm[7] = m[13];  glm[11] = m[14];   glm[15] = m[15];
}

void setLight(GLfloat* a, Vector4& b) {
	a[0] = b.x;
	a[1] = b.y;
	a[2] = b.z;
}

void setReflection(GLfloat* a, Vector3& b) {
	a[0] = b.x;
	a[1] = b.y;
	a[2] = b.z;
}

// Vertex buffers
GLuint VAO, VBO;

// Call back function for window reshape
void ChangeSize(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
	CURRENT_WIDTH = width;
	CURRENT_HEIGHT = height;
	// [TODO] change your aspect ratio
	proj.aspect = (float)(width / 2) / (float)height;
	GLfloat A, B;
	A = (proj.farClip + proj.nearClip) / (proj.nearClip - proj.farClip);
	B = 2 * proj.farClip * proj.nearClip / (proj.nearClip - proj.farClip);
	project_matrix = Matrix4(
		1 / tan(proj.fovy / 120) / proj.aspect, 0.0f, 0.0f, 0.0f,
		0.0f, 1 / tan(proj.fovy / 120), 0.0f, 0.0f,
		0.0f, 0.0f, A, B,
		0.0f, 0.0f, -1.0f, 0.0f
	);
}

// Render function for display rendering
void RenderScene(void) {
	// clear canvas
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	Matrix4 T, R, S, M;
	// [TODO] update translation, rotation and scaling
	T = translate(models[cur_idx].position);
	R = rotate(models[cur_idx].rotation);
	S = scaling(models[cur_idx].scale);
	M = T * R * S;

	Matrix4 MVP;
	GLfloat mvp[16], m[16], v[16], p[16];
	GLfloat ka[3], kd[3], ks[3];
	GLfloat lightPosition[3], lightDirection[3], lightIntensity[3];

	// [TODO] multiply all the matrix
	// row-major ---> column-major
	MVP = project_matrix * view_matrix * M;
	setGLMatrix(mvp, MVP);
	setGLMatrix(m, M);
	setGLMatrix(v, view_matrix);
	setGLMatrix(p, project_matrix);

	if (cur_light_mode == Directional) {
		light_position = Vector4(lightInfo[Directional].position.x, lightInfo[Directional].position.y, lightInfo[Directional].position.z, 0.0f);
		light_direction = Vector4(lightInfo[Directional].direction.x, lightInfo[Directional].direction.y, lightInfo[Directional].direction.z, 0.0f);
		light_intensity = Vector4(lightInfo[Directional].intensity.x, lightInfo[Directional].intensity.y, lightInfo[Directional].intensity.z, 0.0f);
	}
	else if (cur_light_mode == Point) {
		light_position = Vector4(lightInfo[Point].position.x, lightInfo[Point].position.y, lightInfo[Point].position.z, 0.0f);
		light_direction = Vector4(lightInfo[Point].direction.x, lightInfo[Point].direction.y, lightInfo[Point].direction.z, 0.0f);
		light_intensity = Vector4(lightInfo[Point].intensity.x, lightInfo[Point].intensity.y, lightInfo[Point].intensity.z, 0.0f);
	}
	else if (cur_light_mode == Spot) {
		light_position = Vector4(lightInfo[Spot].position.x, lightInfo[Spot].position.y, lightInfo[Spot].position.z, 0.0f);
		light_direction = Vector4(lightInfo[Spot].direction.x, lightInfo[Spot].direction.y, lightInfo[Spot].direction.z, 0.0f);
		light_intensity = Vector4(lightInfo[Spot].intensity.x, lightInfo[Spot].intensity.y, lightInfo[Spot].intensity.z, 0.0f);
	}

	setLight(lightPosition, light_position);
	setLight(lightDirection, light_direction);
	setLight(lightIntensity, light_intensity);

	// use uniform to send mvp to vertex shader
	glUniformMatrix4fv(uniform.iLocMVP, 1, GL_FALSE, mvp);
	glUniformMatrix4fv(uniform.iLocM, 1, GL_FALSE, m);
	glUniformMatrix4fv(uniform.iLocV, 1, GL_FALSE, v);
	glUniformMatrix4fv(uniform.iLocP, 1, GL_FALSE, p);
	glUniform3fv(uniform.iLocPosition, 1, lightPosition);
	glUniform3fv(uniform.iLocDirection, 1, lightDirection);
	glUniform3fv(uniform.iLocIntensity, 1, lightIntensity);
	glUniform1f(uniform.iLocAngle, lightInfo[Spot].angle);
	glUniform1f(uniform.iLocShininess, shininess);
	glUniform1i(uniform.iLocMode, cur_light_mode);
	glUniform1i(uniform.iLocFVMode, 1);
	glViewport(0, 0, CURRENT_WIDTH / 2, CURRENT_HEIGHT);

	for (int i = 0; i < models[cur_idx].shapes.size(); i++)
	{
		setReflection(ka, models[cur_idx].shapes[i].material.Ka);
		setReflection(kd, models[cur_idx].shapes[i].material.Kd);
		setReflection(ks, models[cur_idx].shapes[i].material.Ks);

		glUniform3fv(uniform.iLocKa, 1, ka);
		glUniform3fv(uniform.iLocKd, 1, kd);
		glUniform3fv(uniform.iLocKs, 1, ks);

		// set glViewport and draw twice ... 
		glBindVertexArray(models[cur_idx].shapes[i].vao);
		glDrawArrays(GL_TRIANGLES, 0, models[cur_idx].shapes[i].vertex_count);
	}

	glUniform1i(uniform.iLocFVMode, 0);
	glViewport(CURRENT_WIDTH / 2, 0, CURRENT_WIDTH / 2, CURRENT_HEIGHT);

	for (int i = 0; i < models[cur_idx].shapes.size(); i++)
	{
		// set glViewport and draw twice ... 
		glBindVertexArray(models[cur_idx].shapes[i].vao);
		glDrawArrays(GL_TRIANGLES, 0, models[cur_idx].shapes[i].vertex_count);
	}
}


void KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	// [TODO] Call back function for keyboard
	if (action == GLFW_PRESS) {
		if (key == GLFW_KEY_Z) {	// prev model
			if (cur_idx != 0) {
				cur_idx--;
			}
			else {
				cur_idx = 4;
			}
		}
		else if (key == GLFW_KEY_X) {	// next model
			cur_idx = (cur_idx + 1) % 5;
		}
		else if (key == GLFW_KEY_T) {	// translation
			cur_trans_mode = GeoTranslation;
		}
		else if (key == GLFW_KEY_S) {	// scaling
			cur_trans_mode = GeoScaling;
		}
		else if (key == GLFW_KEY_R) {	// rotation
			cur_trans_mode = GeoRotation;
		}
		else if (key == GLFW_KEY_L) {	// switch light mode
			cur_light_mode = (cur_light_mode + 1) % 3;
		}
		else if (key == GLFW_KEY_K) {	// light edit
			cur_trans_mode = LightEdit;
		}
		else if (key == GLFW_KEY_J) {	// shininess edit
			cur_trans_mode = ShininessEdit;
		}
	}
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	// [TODO] scroll up positive, otherwise it would be negtive
	if (cur_trans_mode == GeoTranslation) {
		models[cur_idx].position.z += yoffset / 6;
	}
	else if (cur_trans_mode == GeoScaling) {
		models[cur_idx].scale.z += yoffset / 6;
	}
	else if (cur_trans_mode == GeoRotation) {
		models[cur_idx].rotation.z += yoffset / 6;
	}
	else if (cur_trans_mode == LightEdit) {
		if (cur_light_mode == Directional) {
			lightInfo[Directional].intensity.x += yoffset / 10;
			lightInfo[Directional].intensity.y += yoffset / 10;
			lightInfo[Directional].intensity.z += yoffset / 10;
		}
		else if (cur_light_mode == Point) {
			lightInfo[Point].intensity.x += yoffset / 10;
			lightInfo[Point].intensity.y += yoffset / 10;
			lightInfo[Point].intensity.z += yoffset / 10;
		}
		else if (cur_light_mode == Spot) {
			lightInfo[Spot].angle += yoffset;
		}
	}
	else if (cur_trans_mode == ShininessEdit) {
		shininess += yoffset;
	}
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	// [TODO] mouse press callback function
	if (action == GLFW_PRESS) {
		double xpos, ypos;
		mouse_pressed = true;
		glfwGetCursorPos(window, &xpos, &ypos);
		starting_press_x = xpos;
		starting_press_y = ypos;
	}
	else if (action == GLFW_RELEASE) {
		mouse_pressed = false;
	}
}

static void cursor_pos_callback(GLFWwindow* window, double xpos, double ypos)
{
	// [TODO] cursor position callback function
	double x_axis, y_axis;
	x_axis = xpos - starting_press_x;
	y_axis = ypos - starting_press_y;
	if (mouse_pressed) {
		if (cur_trans_mode == GeoTranslation) {
			models[cur_idx].position.x += x_axis / WINDOW_WIDTH * 4;
			models[cur_idx].position.y -= y_axis / WINDOW_HEIGHT * 4;
		}
		else if (cur_trans_mode == GeoScaling) {
			models[cur_idx].scale.x -= x_axis / WINDOW_WIDTH;
			models[cur_idx].scale.y -= y_axis / WINDOW_HEIGHT;
		}
		else if (cur_trans_mode == GeoRotation) {
			models[cur_idx].rotation.y -= x_axis / WINDOW_WIDTH * 4;
			models[cur_idx].rotation.x -= y_axis / WINDOW_HEIGHT * 4;
		}
		else if (cur_trans_mode == LightEdit) {
			if (cur_light_mode == Directional) {
				lightInfo[Directional].position.x += x_axis / WINDOW_WIDTH * 4;
				lightInfo[Directional].position.y -= y_axis / WINDOW_WIDTH * 4;
			}
			else if (cur_light_mode == Point) {
				lightInfo[Point].position.x += x_axis / WINDOW_WIDTH * 4;
				lightInfo[Point].position.y -= y_axis / WINDOW_WIDTH * 4;
			}
			else if (cur_light_mode == Spot) {
				lightInfo[Spot].position.x += x_axis / WINDOW_WIDTH * 4;
				lightInfo[Spot].position.y -= y_axis / WINDOW_WIDTH * 4;
			}
		}
		starting_press_x = xpos;
		starting_press_y = ypos;
	}
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
	glAttachShader(p, f);
	glAttachShader(p, v);

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

	uniform.iLocMVP = glGetUniformLocation(p, "matrix.mvp");
	uniform.iLocM = glGetUniformLocation(p, "matrix.m");
	uniform.iLocV = glGetUniformLocation(p, "matrix.v");
	uniform.iLocP = glGetUniformLocation(p, "matrix.p");
	uniform.iLocKa = glGetUniformLocation(p, "reflection.Ka");
	uniform.iLocKd = glGetUniformLocation(p, "reflection.Kd");
	uniform.iLocKs = glGetUniformLocation(p, "reflection.Ks");
	uniform.iLocPosition = glGetUniformLocation(p, "lighting.position");
	uniform.iLocDirection = glGetUniformLocation(p, "lighting.direction");
	uniform.iLocIntensity = glGetUniformLocation(p, "lighting.diffuse_intensity");
	uniform.iLocAngle = glGetUniformLocation(p, "lighting.angle");
	uniform.iLocMode = glGetUniformLocation(p, "cur_light_mode");
	uniform.iLocFVMode = glGetUniformLocation(p, "fragment_mode");
	uniform.iLocShininess = glGetUniformLocation(p, "shininess");

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
	// [TODO] Setup some parameters if you need
	proj.left = -1;
	proj.right = 1;
	proj.top = 1;
	proj.bottom = -1;
	proj.nearClip = 0.001;
	proj.farClip = 100.0;
	proj.fovy = 80;
	proj.aspect = (float)(WINDOW_WIDTH / 2) / (float)WINDOW_HEIGHT;

	main_camera.position = Vector3(0.0f, 0.0f, 2.0f);
	main_camera.center = Vector3(0.0f, 0.0f, 0.0f);
	main_camera.up_vector = Vector3(0.0f, 1.0f, 0.0f);

	lightInfo[Directional].position = Vector3(1.0f, 1.0f, 1.0f);
	lightInfo[Directional].direction = Vector3(0.0f, 0.0f, 0.0f);
	lightInfo[Directional].intensity = Vector3(1.0f, 1.0f, 1.0f);

	lightInfo[Point].position = Vector3(0.0f, 2.0f, 1.0f);
	lightInfo[Point].direction = Vector3(0.0f, 0.0f, 0.0f);
	lightInfo[Point].intensity = Vector3(1.0f, 1.0f, 1.0f);

	lightInfo[Spot].position = Vector3(0.0f, 0.0f, 2.0f);
	lightInfo[Spot].direction = Vector3(0.0f, 0.0f, -1.0f);
	lightInfo[Spot].intensity = Vector3(1.0f, 1.0f, 1.0f);
	lightInfo[Spot].angle = 30;

	setViewingMatrix();
	setPerspective();	//set default projection matrix as perspective matrix
}

void setupRC()
{
	// setup shaders
	setShaders();
	initParameter();

	// OpenGL States and Values
	glClearColor(0.2, 0.2, 0.2, 1.0);
	vector<string> model_list{ "../NormalModels/bunny5KN.obj", "../NormalModels/dragon10KN.obj", "../NormalModels/lucy25KN.obj", "../NormalModels/teapot4KN.obj", "../NormalModels/dolphinN.obj" };
	// [TODO] Load five model at here
	int i = 0;
	while (i < 5) {
		LoadModels(model_list[i]);
		i++;
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
	GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "109000168 HW2", NULL, NULL);
	if (window == NULL)
	{
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);


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
