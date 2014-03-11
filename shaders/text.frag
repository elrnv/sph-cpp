#version 330

uniform sampler2D tex;
uniform vec4 color;

in  vec2 uv;
out vec4 frag_color;

void main()
{
  frag_color = vec4(color.rgb, texture(tex, uv).a);
}
