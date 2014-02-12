#version 330 core

uniform mat4 mvpMtx;
uniform mat4 mvMtx;
uniform mat3 nmlMtx;

in vec3 pos;
in vec3 nml;

// pass vertices to the geometry shader
out fData
{
  vec3 pos;
  vec3 nml;
  vec4 col;
} frag;

void main() 
{
  // populate vertices for geometry shader
  frag.pos = pos; // needed in fragment shader
  frag.nml = normalize(nmlMtx * nml);
  frag.col = vec4(1.0, 1.0, 0.0, 1.0);
  gl_Position = mvpMtx * vec4(pos, 1.0);

}
