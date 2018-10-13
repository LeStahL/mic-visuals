/* Chaos Blocks
 * Copyright (C) 2018 Alexander Kraus <nr4@z10.info>
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
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
 
#version 130

uniform float iTime;
uniform vec2 iResolution;
uniform float iScale;
uniform float iNBeats;
uniform float iHighScale;

const float pi = acos(-1.);
const vec3 c = vec3(1.,0.,-1.);

// Hash function
float rand(vec2 x)
{
    return fract(sin(dot(x-1. ,vec2(12.9898,78.233)))*43758.5453);
}

vec3 rand3(vec3 x)
{
    return vec3(rand(x.x*c.xx),rand(x.y*c.xx),rand(x.z*c.xx));
}

float valuenoise(vec2 x)
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

vec2 add(vec2 sda, vec2 sdb)
{
    return mix(sda, sdb, step(sdb.x, sda.x));
}


vec2 sub(vec2 sda, vec2 sdb)
{
    return mix(-sda, sdb, step(sda.x, sdb.x));
}

// Distance to line segment
float linesegment(vec3 x, vec3 p0, vec3 p1)
{
    vec3 d = p1-p0;
    float t = clamp(dot(x-p0,d)/dot(d,d),0.,1.);
    return length(x-mix(p0,p1,t));
}

// Stroke
float stroke(float sdf, float w)
{
    return abs(sdf)-w;
}

/* compute voronoi distance and closest point.
 * x: coordinate
 * return value: vec3(distance, coordinate of control point)
 */
vec4 vor(vec3 x)
{
    vec3 y = floor(x);
   	float ret = 10.;
    
    //find closest control point. ("In which cell am I?")
    vec3 pf=c.yyy, p;
    float df=100., d;
    
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
            for(int k=-1; k<=1; k+=1)
            {
                p = y + vec3(float(i), float(j), float(k));
                p += rand3(p);

                d = length(x-p);
				
                if(d < df)
                {
                    df = d;
                    pf = p;
                }
            }
    
    //compute voronoi distance: minimum distance to any edge
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
            for(int k=-1; k<=1; k+=1)
            {
                p = y + vec3(float(i), float(j), float(k));
                p += rand3(p);

                vec3 o = p - pf;
                d = abs(.5-dot(x-pf, o)/length(o));
                ret = min(ret, d);
            }
    return vec4(ret, pf);
}

mat3 rot(vec3 p)
{
    return mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))
        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))
        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);
}

vec3 ind;
vec2 scene(vec3 x)
{
    x = rot(.05*vec3(1.,2.,3.)*iTime+iNBeats)*x;
    vec4 v = vor(.5*x-(.1+.05*iScale)*valuenoise(x.xy-2.-1.*iTime)), w  = vor(2.*x-(.2+.1*iScale)*valuenoise(x.xy-2.-1.*iTime));
    ind = v.gba+.1*w.gba;
    float d = max(-stroke(.4*w.x,5.e-3+1.e-3*iScale), stroke(v.x,1.e-1));
    //d = max(-length(x)+1., d);
    return vec2(abs(d)-.001, 1.);
}

vec3 normal(vec3 x)
{
    float dx = 5.e-3;
//    float dx = 5.e-2*clamp(iScale,1.e-1,1.);
    float s = scene(x).x;
    return normalize(vec3(scene(x+dx*c.xyy).x-s, scene(x+dx*c.yxy).x-s, scene(x+dx*c.yyx).x-s));
}

vec3 stdcolor(vec2 x)
{
	return 0.5 + 0.5*cos(iTime+x.xyx+vec3(0,2,4));
}

// Distance to circle
float circle(vec2 x, float r)
{
    return length(x)-r;
}

// Distance to circle segment
float circlesegment(vec2 x, float r, float p0, float p1)
{
    float p = atan(x.y, x.x);
    p = clamp(p, p0, p1);
    return length(x-r*vec2(cos(p), sin(p)));
}

// Distance to line segment
float linesegment(vec2 x, vec2 p0, vec2 p1)
{
    vec2 d = p1-p0;
    float t = clamp(dot(x-p0,d)/dot(d,d),0.,1.);
    return length(x-mix(p0,p1,t));
}

// Distance to 210 logo
float logo(vec2 x, float r)
{
    return min(
        min(circle(x+r*c.zy, r), linesegment(x,r*c.yz, r*c.yx)),
        circlesegment(x+r*c.xy, r, -.5*pi, .5*pi)
    );
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord/iResolution.yy-.5, s;
	vec3 col = c.yyy, ro = c.yyx, r = c.xyy, u = c.yxy, t = c.yxy+ uv.x*r + uv.y*u, x, dir = normalize(t-ro);
    float d = 1.;
    
    for(int i=0; i<80; ++i)
    {
        x = ro + d*dir;
        s = scene(x);
        
        if(s.x < 1.e-4) break;
        if(i == 79)
        {
            col = c.xxx*smoothstep(1.5/iResolution.y, -1.5/iResolution.y, stroke(logo(uv-vec2(-.45,.45),.02),.005));
            fragColor = vec4(col, 1.);
            return;
        }
        
        d += s.x;
    }
    
    vec3 n = normal(x), l = x+c.xxx, re = normalize(reflect(-l,n)), v = normalize(x-ro);
    vec3 c1 = stdcolor(uv+.5*ind.x+iNBeats), 
            c2 = stdcolor(uv+.5*ind.y+iNBeats), 
            c3 = stdcolor(uv+.5*ind.z+iNBeats);
    float rev = abs(dot(re,v)), ln = abs(dot(l,n));
    if(s.y == 1.)
    {
        col = .1*c1*vec3(1.,.3,.3) + .2*c1*vec3(1.,.3,.3)*ln + vec3(1.,1.,.1)*pow(rev,2.*(2.-1.5*clamp(iScale,0.,1.))) + 2.*c1*pow(rev, 8.)+3.*c1*pow(rev, 16.);
        col = abs(col);
    }
    //portability
    col = clamp(.33*col, 0., 1.);
//     col = mix(col, c1*vec3(1.,.3,.3), tanh(1.e-2*length(x-ro)));
    //210 logo
    col = mix(clamp(col,c.yyy,c.xxx), c.xxx, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, stroke(logo(uv-vec2(-.45,.45),.02),.005)));
    //trendy display lines
    col += vec3(0., 0.05, 0.1)*sin(uv.y*1050.+ 5.*iTime);
    
    fragColor = vec4(col,1.0);
}


void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
