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

// pass vertices to the geometry shader
out fData
{
  vec4 pos;
  vec3 nml;
  vec4 col;
} frag;

void main() 
{
  // populate vertices for geometry shader
  frag.pos = ModelMtx * vec4(pos, 1.0f); // needed in fragment shader
  frag.nml = normalize((NormalMtx * vec4(nml, 0.0f)).xyz);
  frag.col = vec4(1.0f, 1.0f, 0.0f, 1.0f);
  gl_Position = MVPMtx * vec4(pos, 1.0f);
}
