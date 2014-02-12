#version 330

in fData
{
  vec4 pos;
  vec3 nml;
  vec4 col;
} frag;

out vec4 frag_color;

void main()
{
   frag_color = frag.col;
}

