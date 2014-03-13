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
  vec4 col;
//  float halo_ratio;
//  float dimming;
} frag;

out vec4 frag_color;

void main()
{
  vec3 acc = vec3(0.0f, 0.0f, 0.0f); // light accumulator

  vec3 E = normalize((eye - frag.pos).xyz);
  vec2 rel = gl_PointCoord - vec2(0.5);

  vec3 nml = normalize((ViewInvMtx * vec4(rel.x, -rel.y, sqrt(1.0f - dot(rel,rel)), 0.0f)).xyz);

  for (int i = 0; i < 2; ++i)
  {
    vec3 L = normalize((lights[i].pos - frag.pos).xyz);
    vec3 R = normalize(reflect(-L, nml));

    vec3 diff = diffuse.xyz * max(dot(nml, L), 0.0f);
    vec3 spec = vec3(0.f, 0.f, 0.f);
    if (options.x > 0)
       spec = specular.xyz * pow( max(dot(R,E), 0.0f), options.x);
    acc += (diff + spec) * lights[i].col.xyz;
  }

//  acc *= frag.dimming;

 if (length(rel) < 0.5)
    frag_color = vec4(ambient.xyz + acc, options.y);
 // if (length(rel) < 0.5*frag.halo_ratio)
 // {
 //   frag_color = vec4(ambient.xyz + acc, options.y);
 // }
 // else if (length(rel) < 0.5 && length(rel) > 0.49)
 // {
 //   frag_color = vec4(0.5*ambient.xyz + acc, 0.5*options.y);
 // }
  else
    discard;
}
