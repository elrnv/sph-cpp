#version 330 core

uniform mat3 mvInvMtx;

uniform vec4 ambientMat;
uniform vec4 diffuseMat;
uniform vec4 specMat;
uniform float specPow;

struct Light
{
  vec3 pos;
  vec4 col;
  vec3 falloff;
  bool camera;
};

uniform Light lights[8];

in fData
{
  vec3 pos;
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
    L = vec3(mvInvMtx * vec4(lights[0].pos, 1.0)) - frag.pos;
  }
  else
  {
    L = lights[0].pos - frag.pos;
  }

  float d2 = dot(L,L);
  vec3 attenvec = vec3(1, sqrt(d), d2);

  float atten = 1 / dot(attenvec, light[0].falloff);

  vec3 E = normalize(-frag.pos);
  vec3 R = normalize(reflect(-L,frag.nml));

  ambient = ambientMat;
  diffuse = diffuseMat * max(dot(frag.nml, normalize(L)), 0.0);
  spec = specMat * pow(max(dot(R,E),0.0),specPow);

  frag_color = ambient + atten * (diffuse + spec);
}
