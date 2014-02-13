#version 330 core

layout (triangles) in;
layout (line_strip, max_vertices=6) out;

uniform mat4 mvpMtx;
uniform mat4 mvMtx;
uniform mat3 nmlMtx;

uniform float normal_length = 0.5;

in vData
{
  vec4 pos;
  vec3 nml;
  vec4 col;
} vtx[];

out fData 
{
  vec4 pos;
  vec3 nml;
  vec4 col;
} frag;

void main()
{
  // wireframe
  int i;
  for (i=0; i < gl_in.length(); ++i)
  {
    frag.pos = vtx[i].pos;
    frag.nml = vtx[i].nml;
    frag.col = vtx[i].col;

    gl_Position = mvpMtx * frag.pos;
    EmitVertex();
  }
  frag.pos = vtx[0].pos;
  frag.nml = vtx[0].nml;
  frag.col = vtx[0].col;

  gl_Position = mvpMtx * frag.pos;
  EmitVertex();

  EndPrimitive();

  // Normals 
  frag.pos = vtx[0].pos;
  frag.nml = vtx[0].nml;
  frag.col = vtx[0].col;

  gl_Position = mvpMtx * frag.pos;
  EmitVertex();

  gl_Position = mvpMtx * vec4(frag.pos.xyz + frag.nml * normal_length, 1.0);
  EmitVertex();

  EndPrimitive();
}
