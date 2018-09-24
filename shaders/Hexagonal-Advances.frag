/*
 * Hexagonal Advances
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
const vec3 c = vec3(1.,0.,-1.);
vec2 ind;

// hash function
float rand(vec2 a0)
{
    return fract(sin(dot(a0.xy ,vec2(12.9898,78.233)))*43758.5453);
}

// value noise
float vnoise(vec2 x)
{
    vec2 y = floor(x);
    x = fract(x);
    float r00 = -1.+2.*rand(y),
        r10 = -1.+2.*rand(y+c.xy),
        r01 = -1.+2.*rand(y+c.yx),
        r11 = -1.+2.*rand(y+c.xx);
    return mix(
        mix(r00, r10, x.x),
        mix(r01, r11, x.x),
        x.y
    );
}

float mfvnoise(vec2 x, float f0, float f1, float phi)
{
    float sum = 0.;
    float a = 1.;
    
    for(float f = f0; f<f1; f = f*2.)
    {
        sum = a*vnoise(f*x) + sum;
        a = a*phi;
    }
    
    return sum;
}

// compute distance to regular star
float dstar(vec2 x, float N, vec2 R)
{
    float d = pi/N,
        p0 = acos(x.x/length(x)),
        p = mod(p0, d),
        i = mod(round((p-p0)/d),2.);
    x = length(x)*vec2(cos(p),sin(p));
    vec2 a = mix(R,R.yx,i),
    	p1 = a.x*c.xy,
        ff = a.y*vec2(cos(d),sin(d))-p1;
   	ff = ff.yx*c.zx;
    return dot(x-p1,ff)/length(ff);
}

// compute distance to regular polygon
float dpoly_min(vec2 x, float N, float R)
{
    float d = 2.*pi/N,
        t = mod(acos(x.x/length(x)), d)-.5*d;
    return R-length(x)*cos(t)/cos(.5*d);
}

float zextrude(float z, float d2d, float h)
{
    vec2 d = abs(vec2(min(d2d, 0.),z))-h*c.yx;
    return min(max(d.x,d.y),0.)+length(max(d,0.));
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

vec2 scene(vec3 x)
{
    x.x += .1*vnoise(x.xy+iTime);
    float dr = .03, d, dh, h;
    vec3 y = vec3(mod(x.xy, dr)-.5*dr, x.z), index = x-y;
    
    vec2 z = (index.xy-c.yx)*(.5);
    float sd = min(cr(z-.125*c.xy, .125, .04), cs(z+.125*c.xy, .125, .04, -pi/2., pi/2.));
    sd = min(sd, b(z, -.125*c.yx, .125*c.yx, .04));
//     sd = abs(sd)-.04;
    
    dh = .05*mfvnoise(index.xy-.5*iTime, 8., 100., .45)+.1*iScale*rand(index.xy);
    h = .1+dh+(.1+iScale*.1)*step(sd,0.);
    d = zextrude(y.z, dstar(rot(y.xy, 5.*iTime), 5.+5.*round(rand(3.*index.xy)), vec2(.2+.1*rand(2.*index.xy),.4+.1*rand(index.xy))*dr), h)-.001;
//     float st = rand(index.xy);
//     if(st > 1.5)
//         d = zextrude(y.z, dpoly_min(y.xy, 6., .5*dr), h)-.001;
//     else
//         
    ind = index.xy;
    
    float guard = -length(max(abs(y)-vec3(.5*dr*c.xx, 20.),0.));
        guard = abs(guard)+dr*.1;
        d = min(d, guard);
        
    vec2 sdf = vec2(d, 1.+step(sd, 0.)), sda = vec2(x.z, 1.);
    
    return mix(sdf, sda, step(sda.x, sdf.x));
}

const float dx = 1.e-3;
vec3 normal(vec3 x)
{
    float s = scene(x).x;
    return normalize(vec3(scene(x+dx*c.xyy).x-s, scene(x+dx*c.yxy).x-s, scene(x+dx*c.yyx).x-s));
}

vec3 background(vec2 x)
{
    return c.yyy;
}

vec3 synthcol(float scale, float phase)
{
    vec3 c2 = .5*vec3(rand(phase*c.xx), rand(phase*c.xx+1.), rand(phase*c.xx+2.))+.5;
    mat3 r1 = rot((5.e-1*phase)*vec3(1.1,1.3,1.5));
    float sc = rand(phase*c.xx);
    if(abs(sc) < .2)
        sc = sign(sc)*.2;
    return 
        (
            sc*1.1*mix
            (
                mix(-(cross(c2, r1*c2)),c.yyy, .5*scale),
                mix(c.yyy, -(r1*c2), .5*scale), 
                scale
            )
        );
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord/iResolution.yy-.5, s;
    vec3 or = rot(c.yyx*iTime)*(c.yyx-.2*c.yxy)+.5*c.yxy, ta = c.yxy, r = c.xyy, u = cross(r, normalize(ta-or)),
        rt = ta+uv.x*r+uv.y*u, rd = normalize(rt-or), x, col;
    float d = 0.;
    
    for(int i=0; i<700; ++i)
    {
        x = or + d * rd;
        s = scene(x);
        
        if(s.x < 1.e-4) break;
        if(i==699)
        {
            fragColor = vec4(background(uv), 1.);
            return;
        }
        
        d += s.x;
    }
    
    vec3 n = normal(x), l = 2.*c.yyx, re = normalize(reflect(-l, n)), v = normalize(or-x);
    
    if(s.y == 1.)
    {
        vec3 c1 = abs(synthcol(1.5+.5*vnoise(2.*(x.xy)), iNBeats+1.1*ind.x*ind.y)), c2 = .1*abs(synthcol(1.5+.5*vnoise(2.*(x.xy)), iNBeats+ind.x*ind.y));
        col = .3*c1+.3*c1*abs(dot(l, n))+c2*dot(re,v);
    }
    else if(s.y == 2.)
    {
        l = c.xxx;
        vec3 c1 = abs(synthcol(1.5+.5*vnoise(2.*(x.xy)), iNBeats+1.1*ind.x*ind.y+20.)), c2 = abs(synthcol(1.5+.5*vnoise(2.*(x.xy)), iNBeats+ind.x*ind.y+20.));
        c1 += .1*mfvnoise(x.xz, 1., 100., .45);
        col = .3*c1+.3*c1*abs(dot(l, n))+c2*dot(re,v);
    }
    
    fragColor = vec4(col, 1.);
}


void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
