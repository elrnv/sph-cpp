#version 330 core

layout(std140) uniform Globals
{
  mat4 MVPMtx;
  mat4 VPMtx;
  mat4 ModelMtx;
  mat4 NormalMtx;
  mat4 ViewInvMtx;
  vec4 eye;
};

in vec3 pos;

void main() 
{
  gl_Position = MVPMtx * vec4(pos, 1.0);
}
