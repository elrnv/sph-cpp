#version 330
in vec3 vertexPosition;
uniform vec3 vertexColor;
uniform mat4 P, V, M;
out vec4 col;
void main() {
	 col = vec4(vertexColor, 1.0);
	 gl_Position = P * V * M * vec4(vertexPosition, 1.0);
}
