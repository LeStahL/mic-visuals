/*
 * Eyeball Explosion
 * 
 * Copyright (C) 2018  Alexander Kraus <nr4@z10.info>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#version 130

uniform float iTime;
uniform vec2 iResolution;
uniform float iScale;
uniform float iNBeats;
uniform float iHighScale;
 
const float pi = acos(-1.);
const vec2 c = vec2(1.,0.);

float rand(vec2 a0)
{
    return fract(sin(dot(a0.xy ,vec2(12.9898,78.233)))*43758.5453);
}

float smoothstep_noise(float x)
{
    float r1 = -1.+2.*rand(floor(x)*c.xx), r2 = -1.+2.*rand(ceil(x)*c.xx);
    return mix(r1, r2, smoothstep(.25, .75, fract(x)));
}

float mfsmoothstep_noise(float x, float f0, float f1, float phi)
{
    float sum = 0.;
    float a = 1.;
    
    for(float f = f0; f<f1; f = f*2.)
    {
        sum = a*smoothstep_noise(f*x) + sum;
        a = a*phi;
    }
    
    return sum;
}

vec2 rot(vec2 x, float p)
{
    return mat2(cos(p), sin(p), -sin(p), cos(p))*x;
}

mat3 rot(vec3 p)
{
    vec3 cp = cos(p), sp = sin(p);
    mat3 m = mat3(cp.y*cp.x, cp.x*sp.z+cp.z*sp.x*sp.y, sp.x*sp.z-cp.x*cp.z*sp.y, 
           -cp.y*sp.z, cp.x*cp.z-sp.x*sp.y*sp.z, cp.z*sp.x+cp.x*sp.y*sp.z, 
           sp.y, -cp.y*sp.x, cp.x*cp.y);
    return m;
}

float rect(vec2 x, vec2 b)
{
    return length(max(abs(x)-b,0.));
}

vec3 synthcol(float scale, float phase)
{
    vec3 c2 = .5*vec3(rand(phase*c.xx), rand(phase*c.xx+1.), rand(phase*c.xx+2.))+.5;
    mat3 r1 = rot((5.e-1*phase)*vec3(1.1,1.3,1.5));
    return 
        (
            .5*rand(phase*c.xx)*1.1*mix
            (
                mix(-(cross(c2, r1*c2)),c.yyy, .5*scale),
                mix(c.yyy, -(r1*c2), .5*scale), 
                scale
            )
        );
}

float cr(vec2 x, float r, float w)
{
    return abs(length(x)-r)-w;
}

float cs(vec2 x, float r0, float w, float p0, float p1)
{
    float r = length(x), p = acos(x.x/r)*step(0.,x.y)-acos(x.x/r)*step(x.y,0.);
    p = clamp(p, p0, p1);
    vec2 y = r0*vec2(cos(p), sin(p));
    return length(x-y)-w;
}

float b(vec2 x, vec2 a, vec2 b, float w)
{
    vec2 d = b-a;
    return length(x-mix(a, b, clamp(dot(x-a, d)/dot(d,d), 0., 1.)))-w;
}

mat3 R(vec3 t)
{
    vec3 ct = cos(t), st = sin(t);
    return mat3(c.xyyy, ct.x, st.x, 0., -st.x, ct.x)
        *mat3(ct.y, 0., -st.y, c.yxy, st.y, 0., ct.y)
        *mat3(ct.z, st.z, 0., -st.z, ct.z, c.yyyx);
}

vec4 add2(vec4 sdf, vec4 sda)
{
    return vec4(
        min(sdf.x, sda.x), 
        mix(sda.gba, sdf.gba, smoothstep(-1.5/iResolution.y, 1.5/iResolution.y, sda.x))
    );
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord/iResolution.yy-.5-.33*c.xy;
    uv = rot(uv, .25*iTime);
    
    vec4 sdf = vec4(0., c.yyy);
    
    float N = 128.,
        r_inner = .2;
    
    for(float i=10.; i>=0.; i-=1.)
    {
        uv = rot(uv, 1.1);
        vec2 p = vec2(length(uv), atan(uv.y/uv.x)-float(i)*.1*iTime),
        q = vec2(p.x-.05*float(i), mod(p.y, 2.*pi/N)-pi/N),
        q0 = vec2(r_inner, q.y);
    
    	float index = (p-q).y;
        
        r_inner = .2+float(i)*.005+ .0005*12.;
        
        float dr = .1*mfsmoothstep_noise(index-iTime-4.*float(i), 1., 100., .45) + .05*rand(index*c.xx+.2*c.yx)+.05*iScale,
            len = abs(.005*float(i)+dr),
            width = abs(.015+.005*rand(index*c.xx+.4));

        vec4 sda = vec4(rect(q-r_inner*c.xy, len*c.xy+width*c.yx), synthcol(((q.x-r_inner)/.05+(q.y/2./pi)),iNBeats));
        sdf = add2(sdf, sda);
    }
    
    //210
    vec2 x = uv*(2.+iScale);
    float sd = min(cr(x-.125*c.xy, .125, .04), cs(x+.125*c.xy, .125, .04, -pi/2., pi/2.));
    sd = min(sd, b(x, -.125*c.yx, .125*c.yx, .04));
    vec4 sda = vec4(sd, 2.*synthcol(1.,iNBeats));
    sdf = add2(sdf, sda);
	vec4 sdb = vec4(abs(sd-.01)-.005, 2.*synthcol(0.,iNBeats));
    sdf = add2(sdf, sdb);

    vec3 col = sdf.gba*smoothstep(1.5/iResolution.y, -1.5/iResolution.y, sdf.x);

    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
