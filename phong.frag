#version 330 core

uniform mat4 mvInvMtx;
uniform mat3 nmlMtx;
uniform mat4 mvMtx;

uniform vec4 ambientMat;
uniform vec4 diffuseMat;
uniform vec4 specMat;
uniform float specPow;

struct Light
{
  vec4 pos;
  vec4 col;
  bool camera;
};

uniform Light lights[8];

in fData
{
  vec4 pos;
  vec3 nml;
  vec4 col;
} frag;

out vec4 frag_color;

void main()
{
  vec4 diffuse;
  vec4 spec;
  vec4 ambient;

  vec3 L; // vector from vtx to light source
  if (lights[0].camera)
  {
    L = normalize(vec3(mvInvMtx * lights[0].pos - frag.pos));
  }
  else
  {
    L = normalize(vec3(lights[0].pos - frag.pos));
  }

  vec3 E = normalize(-vec3(mvMtx*frag.pos));
  vec3 R = normalize(reflect(-L, frag.nml));

  ambient = ambientMat;
  diffuse = diffuseMat * max(dot(frag.nml, L), 0.0f);
  spec = specMat * pow( max(dot(R,E), 0.0f), specPow );

  frag_color = ambient + (diffuse + spec);// * lights[0].col;
}
