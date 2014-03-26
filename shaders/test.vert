#version 330

in vec2 vtx;

void main()
{
  gl_Position = vec4((vec2(10, 10) + vtx.xy)*vec2(2.0/600.0, -2.0/800.0)+ vec2(-1.0f, 1.0f), 0.0f, 1.0f);
}
