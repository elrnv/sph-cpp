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

// pass vertices to the geometry shader
out vData
{
  vec4 col;
} vtx;

void main() 
{
  // populate vertices for geometry shader
  vtx.col = vec4(1.0f, 1.0f, 0.0f, 1.0f);
  gl_Position = ModelMtx * vec4(pos, 1.0f);
}
