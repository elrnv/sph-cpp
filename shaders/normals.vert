#version 330 core

layout(std140) uniform Globals
{
  mat4 MVPMtx;
  mat4 VPMtx;
  mat4 ModelMtx;
  mat4 NormalMtx;
  mat4 ViewInvMtx;
  vec4 eye;

  vec4 ambient;
  vec4 diffuse;
  vec4 specular;
  vec4 options;
};

in vec3 pos;
in vec3 nml;

out vec3 vtxnml;

void main() 
{
  vtxnml = normalize((NormalMtx * vec4(nml, 0.0f)).xyz);
  gl_Position = ModelMtx * vec4(pos, 1.0);
}
