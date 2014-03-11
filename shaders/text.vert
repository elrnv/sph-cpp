#version 330

uniform float xpos;
uniform float ypos;
uniform vec2 scale;

in vec4 vtx;
out vec2 uv;

void main()
{
  uv = vtx.zw;
  gl_Position = vec4((vec2(xpos, ypos) + vtx.xy)*scale + vec2(-1.0f, 1.0f), 0.0f, 1.0f);
}
