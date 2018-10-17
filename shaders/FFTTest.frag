/*
 * 210 Logo
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
uniform sampler2D iFFT;
uniform float iFFTWidth;

/*************************************/
/** Constants and utility functions **/
/*************************************/

// Constants
const float pi = acos(-1.);
const vec3 c = vec3(1.,0.,-1.);

// Image rescaling function
float rescale(float value, vec2 old, vec2 new)
{
    return mix(new.x, new.y, (value-old.x)/(old.y-old.x));
}

/**********************/
/** Color generation **/
/**********************/

// Standard shadertoy color
vec3 stdcolor(vec2 x)
{
	return 0.5 + 0.5*cos(iTime+x.xyx+vec3(0,2,4));
}

/***********/
/** Noise **/
/***********/

// Hash function
float rand(vec2 x)
{
    return fract(sin(dot(x-1. ,vec2(12.9898,78.233)))*43758.5453);
}

// 2D value noise
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

// Multi-frequency value noise
float mfvaluenoise(vec2 x, float f0, float f1, float phi)
{
    float sum = 0.;
    float a = 1.2;
    
    for(float f = f0; f<f1; f = f*2.)
    {
        sum = a*valuenoise(f*x) + sum;
        a = a*phi;
    }
    
    return sum;
}

/****************/
/** 2D Objects **/
/****************/

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

/**********************/
/** Scene management **/
/**********************/

// Add objects to scene with proper antialiasing
vec4 add(vec4 sdf, vec4 sda)
{
    return vec4(
        min(sdf.x, sda.x), 
        mix(sda.gba, sdf.gba, smoothstep(-1.5/iResolution.y, 1.5/iResolution.y, sda.x))
    );
}

/******************/
/** 2D operators **/
/******************/

// Add width to object
float stroke(float sdf, float w)
{
    return abs(sdf)-w;
}

// Add thickness effect to object
vec4 thick(vec2 x, vec4 sdf, vec2 n)
{
    for(int i=1; i<6; ++i)
		sdf = add(vec4(stroke(sdf.x*n.x*n.y*2.*valuenoise((3.+4.*iScale)*x-2.-1.*iTime-1.2), .01), 3.e-3/abs(sdf.x+.2*valuenoise(x-2.-1.*iTime))*stdcolor(x+c.xx*.3*float(i))), sdf); 
    return sdf;
}

// Draw Geometry
vec4 geometry(vec2 x)
{
    vec4 sdf = vec4(stroke(stroke(logo(x, .2), .06),.01), 2.5*stdcolor(x*1.7));
    //for(int i=0; i<10; ++i)
    //    sdf = add(sdf, vec4(stroke(circle(x-.5*vec2(valuenoise(x-2.-5.*iTime+2.*rand(float(i+3)*c.xx)), valuenoise(x-2.-5.*iTime+rand(float(i)*c.xx)))-.5*c.xy, .2+valuenoise(x-2.-5.*iTime+rand(float(i)*c.xx))),.01), 2.5*stdcolor(x+float(i)*.1)));
    return sdf;
}

// Normal
const float dx = 1.e-4;
vec2 normal(vec2 x)
{
    float s = geometry(x).x;
    return normalize(vec2(geometry(x+dx*c.xy).x-s, geometry(x+dx*c.yx).x-s));
}

// Add Effects
vec4 scene(vec2 x)
{
    // Preprocessing
    x += .1*vec2(valuenoise(x-5.*iTime), valuenoise(x-2.-5.*iTime));
    
    // Normal scene
    vec4 sdf = geometry(x);
    vec2 n = normal(x);
    
    // Add potential effect
    sdf = thick(x, sdf, n);
    
    // Cover artefacts
    
    return sdf;
}

float value(float i)
{
//     i /= iFFTWidth;
    //convert 1d index to 2d index and map to texture coordinates
    vec2 ixy = vec2(floor(i/iFFTWidth), mod(i, iFFTWidth));
    //rescale to normal values
    return dot(texture(iFFT, ixy.yx/iFFTWidth), vec4(5.960464477539063e-8, 0.0000152587890625, 0.00390625, 1.));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord/iResolution.yy-.5;
    float v = value(iFFTWidth*iFFTWidth*(uv.x+.5));
    vec3 col = step(0., 2.*v-uv.y-.5)*c.xyy;
//     vec3 col = c.xyy * dot(texture(iFFT, uv+.5),vec4(5.960464477539063e-8, 0.0000152587890625, 0.00390625, 1.));
    /*vec4 s;
    //if(iTime < 20.)
    {
        s = scene(uv);
    }
    //else if(iTime < 40.)
    {
        
    }
    vec3 col = s.gba * smoothstep(1.5/iResolution.y, -1.5/iResolution.y, s.x);
    
    //small 210 logo
    col = mix(clamp(col,c.yyy,c.xxx), c.xxx, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, stroke(logo(uv-vec2(-.45,.45),.02),.005)));
    
    col += vec3(0., 0.05, 0.1)*sin(uv.y*1050.+ 5.*iTime);
*/
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
