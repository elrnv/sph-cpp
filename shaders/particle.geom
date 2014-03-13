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

// radius of the halo of a point
uniform float pt_halo;

in vData
{
  vec4 col;
} vtx[];

out fData 
{
  vec4 pos;
  vec4 col;
  float halo_ratio;
  float dimming;
} frag;

void main()
{
  frag.pos = gl_in[0].gl_Position;
  frag.col = vtx[0].col;
  frag.halo_ratio = pt_radius / pt_halo;
  gl_Position = VPMtx * frag.pos;
  gl_PointSize = 2.0 * pt_halo * pt_scale / distance(eye, frag.pos);
  frag.dimming = 2.0 / distance(eye, frag.pos);
  EmitVertex();
  EndPrimitive();
}
