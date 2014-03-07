#version 330 core

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

struct Light
{
  vec4 pos;
  vec4 col;
};

uniform Light lights[2];

in fData
{
  vec4 pos;
  vec3 nml;
  vec4 col;
} frag;

out vec4 frag_color;

void main()
{
  vec3 acc = vec3(0.0f, 0.0f, 0.0f); // light accumulator

  vec3 E = normalize((eye - frag.pos).xyz);

  for (int i = 0; i < 2; ++i)
  {
    vec3 L = normalize((lights[i].pos - frag.pos).xyz);
    vec3 R = normalize(reflect(-L, frag.nml));

    vec3 diff = diffuse.xyz * max(dot(frag.nml, L), 0.0f);
    vec3 spec = vec3(0.f, 0.f, 0.f);
    if (options.x > 0)
       spec = specular.xyz * pow( max(dot(R,E), 0.0f), options.x);

    acc += (diff + spec) * lights[i].col.xyz;
  }

  frag_color = vec4(ambient.xyz + acc, options.y);
}
