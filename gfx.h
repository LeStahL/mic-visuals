/* File generated with Shader Minifier 1.1.5
 * http://www.ctrl-alt-test.fr
 */
#ifndef GFX_H_
# define GFX_H_

const char *gfx_frag =
 "#version 130\n"
 "uniform float iTime;"
 "uniform vec2 iResolution;"
 "uniform float iScale,iNBeats,iHighScale;"
 "const float pi=acos(-1.);"
 "const vec2 c=vec2(1.,0.);"
 "float rand(vec2 v)"
 "{"
   "return fract(sin(dot(v.xy,vec2(12.9898,78.233)))*43758.5);"
 "}"
 "float smoothstep_noise(float v)"
 "{"
   "float y=-1.+2.*rand(floor(v)*c.xx),f=-1.+2.*rand(ceil(v)*c.xx);"
   "return mix(y,f,smoothstep(.25,.75,fract(v)));"
 "}"
 "float mfsmoothstep_noise(float y,float v,float x,float s)"
 "{"
   "float f=0.,a=1.;"
   "for(float r=v;r<x;r=r*2.)"
     "f=a*smoothstep_noise(r*y)+f,a=a*s;"
   "return f;"
 "}"
 "vec2 rot(vec2 y,float v)"
 "{"
   "return mat2(cos(v),sin(v),-sin(v),cos(v))*y;"
 "}"
 "mat3 rot(vec3 v)"
 "{"
   "vec3 f=cos(v),s=sin(v);"
   "mat3 y=mat3(f.y*f.x,f.x*s.z+f.z*s.x*s.y,s.x*s.z-f.x*f.z*s.y,-f.y*s.z,f.x*f.z-s.x*s.y*s.z,f.z*s.x+f.x*s.y*s.z,s.y,-f.y*s.x,f.x*f.y);"
   "return y;"
 "}"
 "float rect(vec2 v,vec2 y)"
 "{"
   "return length(max(abs(v)-y,0.));"
 "}"
 "vec3 synthcol(float y,float v)"
 "{"
   "vec3 s=.5*vec3(rand(v*c.xx),rand(v*c.xx+1.),rand(v*c.xx+2.))+.5;"
   "mat3 f=rot(.5*v*vec3(1.1,1.3,1.5));"
   "return.5*rand(v*c.xx)*1.1*mix(mix(-cross(s,f*s),c.yyy,.5*y),mix(c.yyy,-(f*s),.5*y),y);"
 "}"
 "float cr(vec2 v,float y,float x)"
 "{"
   "return abs(length(v)-y)-x;"
 "}"
 "float cs(vec2 v,float y,float f,float s,float x)"
 "{"
   "float i=length(v),r=acos(v.x/i)*step(0.,v.y)-acos(v.x/i)*step(v.y,0.);"
   "r=clamp(r,s,x);"
   "vec2 m=y*vec2(cos(r),sin(r));"
   "return length(v-m)-f;"
 "}"
 "float b(vec2 v,vec2 y,vec2 f,float x)"
 "{"
   "vec2 s=f-y;"
   "return length(v-mix(y,f,clamp(dot(v-y,s)/dot(s,s),0.,1.)))-x;"
 "}"
 "mat3 R(vec3 v)"
 "{"
   "vec3 f=cos(v),s=sin(v);"
   "return mat3(c.xyyy,f.x,s.x,0.,-s.x,f.x)*mat3(f.y,0.,-s.y,c.yxy,s.y,0.,f.y)*mat3(f.z,s.z,0.,-s.z,f.z,c.yyyx);"
 "}"
 "vec4 add2(vec4 f,vec4 v)"
 "{"
   "return vec4(min(f.x,v.x),mix(v.yzw,f.yzw,smoothstep(-1.5/iResolution.y,1.5/iResolution.y,v.x)));"
 "}"
 "void mainImage(out vec4 y,in vec2 v)"
 "{"
   "vec2 f=v/iResolution.yy-.5-.33*c.xy;"
   "f=rot(f,.25*iTime);"
   "vec4 s=vec4(0.,c.yyy);"
   "float x=128.,r=.2;"
   "for(float m=10.;m>=0.;m-=1.)"
     "{"
       "f=rot(f,1.1);"
       "vec2 a=vec2(length(f),atan(f.y/f.x)-float(m)*.1*iTime),i=vec2(a.x-.05*float(m),mod(a.y,2.*pi/x)-pi/x),p=vec2(r,i.y);"
       "float l=(a-i).y;"
       "r=.2+float(m)*.005+.006;"
       "float e=.1*mfsmoothstep_noise(l-iTime-4.*float(m),1.,100.,.45)+.05*rand(l*c.xx+.2*c.yx)+.05*iScale,z=abs(.005*float(m)+e),g=abs(.015+.005*rand(l*c.xx+.4));"
       "vec4 b=vec4(rect(i-r*c.xy,z*c.xy+g*c.yx),synthcol((i.x-r)/.05+i.y/2./pi,iNBeats));"
       "s=add2(s,b);"
     "}"
   "vec2 m=f*(2.+iScale);"
   "float a=min(cr(m-.125*c.xy,.125,.04),cs(m+.125*c.xy,.125,.04,-pi/2.,pi/2.));"
   "a=min(a,b(m,-.125*c.yx,.125*c.yx,.04));"
   "vec4 i=vec4(a,2.*synthcol(1.,iNBeats));"
   "s=add2(s,i);"
   "vec4 e=vec4(abs(a-.01)-.005,2.*synthcol(0.,iNBeats));"
   "s=add2(s,e);"
   "vec3 l=s.yzw*smoothstep(1.5/iResolution.y,-1.5/iResolution.y,s.x);"
   "y=vec4(l,1.);"
 "}"
 "void main()"
 "{"
   "mainImage(gl_FragColor,gl_FragCoord.xy);"
 "}";

#endif // GFX_H_
