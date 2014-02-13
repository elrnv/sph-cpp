#version 330 core

uniform mat4 mvpMtx;

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
  frag.pos = vec4(pos, 1.0); // needed in fragment shader
  frag.nml = nml;
  frag.col = vec4(1.0, 1.0, 0.0, 1.0);
  gl_Position = mvpMtx * frag.pos;

}
