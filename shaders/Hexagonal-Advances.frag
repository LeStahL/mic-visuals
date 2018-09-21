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

// hash function
float r(vec2 a0)
{
    return fract(sin(dot(a0.xy ,vec2(12.9898,78.233)))*43758.5453);
}

// compute distance to regular polygon
float dpoly_min(vec2 x, float N, float R)
{
    float d = 2.*pi/N,
        t = mod(acos(x.x/length(x)), d)-.5*d;
    return R-length(x)*cos(t)/cos(.5*d);
}

vec2 scene(vec3 x)
{
    return vec2(length(x-c.yxy)-.1, 1.);
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

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord/iResolution.yy-.5, s;
    vec3 or = c.yyx, ta = c.yxy, r = c.xyy, u = cross(r, normalize(ta-or)),
        rt = ta+uv.x*r+uv.y*u, rd = normalize(rt-or), x, col;
    float d = 0.;
    
    for(int i=0; i<100; ++i)
    {
        x = or + d * rd;
        s = scene(x);
        
        if(s.x < 1.e-4) break;
        if(i==99)
        {
            fragColor = vec4(background(uv), 1.);
            return;
        }
        
        d += s.x;
    }
    
    vec3 n = normal(x), l = c.yyx, re = normalize(reflect(-l, n)), v = normalize(or-x);
    
    if(s.y == 1.)
    {
        col = .1*c.xyy+.3*c.xyy*dot(l, n)+c.xxy*dot(re,v);
    }
    
    fragColor = vec4(col, 1.);
}


void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
