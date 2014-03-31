#version 330 core

layout(std140) uniform Globals
{
  mat4 MVPMtx;
  mat4 VPMtx;
  mat4 ModelMtx;
  mat4 NormalMtx;
  mat4 ViewInvMtx;
  vec4 eye;
};

uniform vec3 ambient;
uniform vec3 diffuse;
uniform vec3 specular;
uniform float specpow;
uniform float opacity;

struct Light
{
  vec4 pos;
  vec4 col;
};

uniform Light lights[2];

in fData
{
  vec4 pos;
  vec4 col;
  float halo_ratio;
  float dimming;
} frag;

out vec4 frag_color;

void main()
{
  vec3 acc = vec3(0.0f, 0.0f, 0.0f); // light accumulator

  vec3 E = normalize((eye - frag.pos).xyz);
  vec2 rel_whole = 2.0*(gl_PointCoord - vec2(0.5));
  vec2 rel = 0.9*(1/frag.halo_ratio)*rel_whole;

  float len_rel = length(rel_whole);
  if (len_rel < frag.halo_ratio)
  {
    vec3 nml = normalize((ViewInvMtx * vec4(rel.x, -rel.y,
            sqrt(1.0f - dot(rel,rel)), 0.0f)).xyz);

    for (int i = 0; i < 2; ++i)
    {
      vec3 L = normalize((lights[i].pos - frag.pos).xyz);
      vec3 R = normalize(reflect(-L, nml));

      vec3 diff = diffuse * max(dot(nml, L), 0.0f);
      vec3 spec = vec3(0.f, 0.f, 0.f);
      if (specpow > 0)
         spec = specular * pow( max(dot(R,E), 0.0f), specpow);
      acc += (diff + spec) * lights[i].col.xyz;
    }

    acc *= frag.dimming;

    frag_color = vec4((ambient + acc)*opacity, 1.0);
  }
  else if (len_rel < 1 && len_rel > 0.95)
  {
    frag_color = vec4(0.1*(ambient + lights[0].col.xyz)*opacity, 1.0);
  }
  else
    discard;
}
