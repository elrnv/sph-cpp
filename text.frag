#version 330

uniform sampler2D tex;
uniform vec4 color;
 
in vec2 texcoord;
out vec4 fragColor;

void main()
{
  fragColor = vec4(color.rgb, texture2D(tex, texcoord).a);
}
