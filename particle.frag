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

struct Light
{
  vec4 pos;
  vec4 col;
};

uniform Light lights[2];
uniform vec4 ambientColor;

in fData
{
  vec4 pos;
  vec4 col;
} frag;

uniform vec4 diff;
uniform vec4 spec;
uniform float specpow;

out vec4 frag_color;

void main()
{
  vec4 acc = vec4(0.0f, 0.0f, 0.0f, 0.0f); // light accumulator

  vec3 E = normalize((eye - frag.pos).xyz);
  vec2 rel = gl_PointCoord - vec2(0.5);

  vec3 nml = normalize((ViewInvMtx * vec4(rel.x, -rel.y, sqrt(1.0f - dot(rel,rel)), 0.0f)).xyz);

  for (int i = 0; i < 2; ++i)
  {
    vec3 L = normalize((lights[i].pos - frag.pos).xyz);
    vec3 R = normalize(reflect(-L, nml));

    vec4 diffuse = diff * max(dot(nml, L), 0.0f);
    vec4 specular = spec * pow( max(dot(R,E), 0.0f), specpow);
    acc += (diffuse + specular) * lights[i].col;
  }

  if (length(rel) < 0.5)
  {
    frag_color = ambientColor + acc;
  }
  else
  {
    frag_color = vec4(0.0f,0.0f,0.0f, 0.0f);
  }
}
