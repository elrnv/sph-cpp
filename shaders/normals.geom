#version 330 core

layout (triangles) in;
layout (line_strip, max_vertices=6) out;

layout(std140) uniform Globals
{
  mat4 MVPMtx;
  mat4 VPMtx;
  mat4 ModelMtx;
  mat4 NormalMtx;
  mat4 ViewInvMtx;
  vec4 eye;
};

uniform float normal_length = 0.1;
uniform vec4  wirecolor = vec4(0.2, 0.5, 0.2, 1.0);

in vec3 vtxnml[];

uniform vec3 diffuse;
uniform float opacity;

out vec4 color;

void main()
{
  // wireframe
  for (int i = 0; i < 3; ++i)
  {
    color = vec4(diffuse, opacity);
    gl_Position = VPMtx * gl_in[i].gl_Position;
    EmitVertex();
  }

  color = vec4(diffuse, opacity);

  gl_Position = VPMtx * gl_in[0].gl_Position;
  EmitVertex();

  EndPrimitive();

  // Normals 
  color = vec4(1.0, 0.0, 0.0, 1.0);

  gl_Position = VPMtx * gl_in[0].gl_Position;
  EmitVertex();

  color = vec4(1.0, 0.0, 0.0, 1.0);

  gl_Position = VPMtx * vec4(gl_in[0].gl_Position.xyz + vtxnml[0] * normal_length, 1.0);
  EmitVertex();

  EndPrimitive();
}
