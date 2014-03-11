#version 330 core

layout (points) in;
layout (points, max_vertices=1) out;

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

// canonically this is near plane distance from eye in pixels
uniform float pt_scale;

// radius of point in world space
uniform float pt_radius;

in vData
{
  vec4 col;
} vtx[];

out fData 
{
  vec4 pos;
  vec4 col;
} frag;

void main()
{
  frag.pos = gl_in[0].gl_Position;
  frag.col = vtx[0].col;
  gl_Position = VPMtx * frag.pos;
  gl_PointSize = 2.0 * pt_radius * pt_scale / distance(eye, frag.pos);
  EmitVertex();
  EndPrimitive();
}
