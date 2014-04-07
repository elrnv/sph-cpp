#version 330

out vec4 frag_color;

uniform vec3 diffuse;
uniform float opacity;

void main()
{
   frag_color = vec4(diffuse, opacity);
}
