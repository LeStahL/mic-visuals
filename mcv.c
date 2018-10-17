/* The Team210 mic to visuals tool for DJs
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
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

int _fltused = 0;

#define ABS(x) ((x)<0?(-x):(x))

#define WIN32_LEAN_AND_MEAN
#define VC_EXTRALEAN
#include <windows.h>
#include <mmsystem.h>

#include <GL/gl.h>
#include "glext.h"

#include <complex.h> 
#include "fftw3.h"

#include <math.h>

//TODO: remove
#include <stdio.h>

// Standard library and CRT rewrite for saving executable size
void *memset(void *ptr, int value, size_t num)
{
    for(int i=num-1; i>=0; i--)
        ((unsigned char *)ptr)[i] = value;
    return ptr;
}

size_t strlen(const char *str)
{
    int len = 0;
    while(str[len] != '\0') ++len;
    return len;
}

void *malloc( unsigned int size )
{
    return GlobalAlloc(GMEM_ZEROINIT, size) ;
}

// OpenGL extensions
PFNGLGETPROGRAMIVPROC glGetProgramiv;
PFNGLGETSHADERIVPROC glGetShaderiv;
PFNGLGETSHADERINFOLOGPROC glGetShaderInfoLog;
PFNGLCREATESHADERPROC glCreateShader;
PFNGLCREATEPROGRAMPROC glCreateProgram;
PFNGLSHADERSOURCEPROC glShaderSource;
PFNGLCOMPILESHADERPROC glCompileShader;
PFNGLATTACHSHADERPROC glAttachShader;
PFNGLLINKPROGRAMPROC glLinkProgram;
PFNGLUSEPROGRAMPROC glUseProgram;
PFNGLGETUNIFORMLOCATIONPROC glGetUniformLocation;
PFNGLUNIFORM2FPROC glUniform2f;
PFNGLUNIFORM1FPROC glUniform1f;
PFNGLGENFRAMEBUFFERSPROC glGenFramebuffers;
PFNGLBINDFRAMEBUFFERPROC glBindFramebuffer;
PFNGLFRAMEBUFFERTEXTURE2DPROC glFramebufferTexture2D;
PFNGLNAMEDRENDERBUFFERSTORAGEEXTPROC glNamedRenderbufferStorageEXT;
PFNGLACTIVETEXTUREPROC glActiveTexture;
PFNGLUNIFORM1IPROC glUniform1i;

void debug(int shader_handle)
{
    int compile_status = 0;
    glGetShaderiv(shader_handle, GL_COMPILE_STATUS, &compile_status);
    if(compile_status != GL_TRUE)
    {
        printf("FAILED.\n");
        int len = 12;
        glGetShaderiv(shader_handle, GL_INFO_LOG_LENGTH, &len);
        printf("log length: %d\n", len);
        GLchar CompileLog[1024];
        glGetShaderInfoLog(shader_handle, len, NULL, CompileLog);
        printf("error: %s\n", CompileLog);
    }
    else 
        printf("shader compilation successful.\n");
}

//Shader globals
int w = 1920, h = 1080,
    index=0, nfiles=0,
    *handles,*programs,*time_locations,
    *resolution_locations,*scale_locations,*nbeats_locations,
    *highscale_locations, *fft_texture_locations, *fft_width_locations;
    
//Demo globals
float t_start = 0., scale = 0., max = -1., nbeats = 0., highscale=0.;
int samplerate = 44100;
WAVEHDR headers[2];
HWAVEIN wi;
int cutoff = 8;

//Recording globals
int double_buffered = 0;
int buffer_size = 512;

//FFTW3 globals
#define NFFT 2048
unsigned int fft_texture_handle;
int fft_texture_size, fft_texture_location, fft_texture_width_location;
float values[NFFT];
fftw_complex *in, *out;
fftw_plan p;
float power_spectrum[NFFT];

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    switch(uMsg)
    {
        case WM_KEYDOWN:
            switch(wParam)
            {
                case VK_ESCAPE:
                    ExitProcess(0);
                    break;
                    
                case VK_LEFT:
                    --index;
                    if(index <0)index = nfiles-1;
                    glUseProgram(programs[index]);
                    printf("index %d", index);
                    break;
                    
                case VK_RIGHT:
                    ++index;
                    if(index == nfiles)index = 0;
                    glUseProgram(programs[index]);
                    printf("index %d", index);
                    break;
                
                case VK_CONTROL:
                    scale = 1.;
                    nbeats += 1.;
                    break;
                    
                case VK_DOWN:
                    cutoff = max(cutoff*3/4,3);
                    break;
                    
                case VK_UP:
                    cutoff = min(cutoff*4/3, NFFT-1);
                    break;
            }
            break;
            
        case WM_TIMER:
            for(int i=0; i<double_buffered+1; ++i)
            {
                if(headers[i].dwFlags & WHDR_DONE)
                {
                    // replace last block in values
                    for(int j=0; j<NFFT-buffer_size; ++j)
                        values[j] = values[j+buffer_size];
                    
                    for(int j=0; j<buffer_size; ++j)
//                         values[NFFT-buffer_size+j] = ((float*)headers[i].lpData)[j];
                        values[NFFT-buffer_size+j] = ((float)(*(short *)(headers[i].lpData+2*j))/32767.);
                    
//                     for(int j=0; j<NFFT; ++j)
//                         printf("%d %le\n", j, values[j]);
                    
                    // fourier transform values
                    for(int j=0; j<NFFT; ++j)
                    {
//                         int low = max(0, j-10);
//                         int high = min(NFFT, j+10);
//                         in[j][0] = 0.;
//                         for(int k=low; k<high; ++k)
//                             in[j][0] += values[k] * (.54 - .46*cos(2.*acos(-1.)*(float)k/(float)NFFT));
                        in[j][0] = values[j] * (.54 - .46*cos(2.*acos(-1.)*(float)j/((float)NFFT-1.)));
                        in[j][1] = 0.;
                    }
                    fftw_execute(p);
                    
                    // compute uniform contents
                    scale = 0.;
                    highscale = 0.;
                    float pmax = 0., pmin = 1.e9;
                    for(int j=0; j<NFFT; ++j)
                    {
                        power_spectrum[j] = out[j][0]*out[j][0]+out[j][1]*out[j][1];
                        pmax = max(pmax, power_spectrum[j]);
                        pmin = min(pmin, power_spectrum[j]);
//                         pm += power_spectrum[j];
                    }
                    for(int j=0; j<NFFT; ++j)
                    {
                        power_spectrum[j] -= pmin;
                        power_spectrum[j] = max(power_spectrum[j], 0.);
                        power_spectrum[j] = min(power_spectrum[j], 1.);
                        //power_spectrum[j] /= (pmax-pmin);
//                         printf("%le\n", power_spectrum[j]);
                    }
                    
                    for(int j=0; j<cutoff; ++j)
                    {
                        scale += power_spectrum[j];
                    }

                    for(int j=cutoff; j<NFFT; ++j)
                    {
                        highscale += power_spectrum[j];
                    }
                    
//                     scale/=(float)cutoff;
//                     highscale/=(float)(NFFT-cutoff);
                    
//                     scale *= 1.e2;
//                     highscale *= 1.e2;
//                     scale *= 2e1;
//                     printf("%le %le\n", scale, .5);
//                     highscale *= .5;
                    
                    if(scale > .5)
                        nbeats += 1.;
                    
                    headers[i].dwFlags = 0;
                    headers[i].dwBytesRecorded = 0;

                    waveInPrepareHeader(wi, &headers[i], sizeof(headers[i]));
                    waveInAddBuffer(wi, &headers[i], sizeof(headers[i]));
                    
                    
                }
            }
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, fft_texture_size, fft_texture_size, 0, GL_RGBA, GL_FLOAT, power_spectrum);
            HDC hdc = GetDC(hwnd);
            
            glUniform1i(fft_texture_locations[index], 0);
            glUniform1f(fft_width_locations[index], fft_texture_size);
            
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, fft_texture_handle);
            
            SYSTEMTIME st_now;
            GetSystemTime(&st_now);
            float t_now = (float)st_now.wMinute*60.+(float)st_now.wSecond+(float)st_now.wMilliseconds/1000.;
            
            glUniform1f(time_locations[index], t_now-t_start);
            glUniform2f(resolution_locations[index], w, h);
            glUniform1f(scale_locations[index], scale);
            glUniform1f(nbeats_locations[index], nbeats);
            glUniform1f(highscale_locations[index], highscale);
            
            glBegin(GL_QUADS);
            
            glVertex3f(-1,-1,0);
            glVertex3f(-1,1,0);
            glVertex3f(1,1,0);
            glVertex3f(1,-1,0);
            
            glEnd();

            glFlush();
            
            SwapBuffers(hdc);
            
            break;
            
        default:
            break;
            
    }
    return DefWindowProc(hwnd, uMsg, wParam, lParam);
}

int WINAPI demo(HINSTANCE hInstance, HINSTANCE hPrevInstance, PWSTR pCmdLine, int nCmdShow)
{
    //TODO: remove
    AllocConsole();
    freopen("CONIN$", "r", stdin);
    freopen("CONOUT$", "w", stdout);
    freopen("CONOUT$", "w", stderr);
    
    CHAR WindowClass[]  = "Team210 Demo Window";
    
    WNDCLASSEX wc = { 0 };
    wc.cbSize = sizeof(wc);
    wc.style = CS_OWNDC | CS_VREDRAW | CS_HREDRAW;
    wc.lpfnWndProc = &WindowProc;
    wc.cbClsExtra = 0;
    wc.cbWndExtra = 0;
    wc.hInstance = hInstance;
    wc.hIcon = LoadIcon(NULL, IDI_WINLOGO); 
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.hbrBackground = NULL;
    wc.lpszMenuName = NULL;
    wc.lpszClassName = WindowClass;
    wc.hIconSm = NULL;
    
    RegisterClassEx(&wc);
    
    // Get full screen information
    HMONITOR hmon = MonitorFromWindow(0, MONITOR_DEFAULTTONEAREST);
    MONITORINFO mi = { sizeof(mi) };
    GetMonitorInfo(hmon, &mi);
    
    // Create the window.
    HWND hwnd = CreateWindowEx(
        0,                                                          // Optional window styles.
        WindowClass,                                                // Window class
        ":: NR4^QM/Team210 :: GO - MAKE A DEMO ::",                                 // Window text
        WS_POPUP | WS_VISIBLE,                                      // Window style
        mi.rcMonitor.left,
        mi.rcMonitor.top,
        mi.rcMonitor.right - mi.rcMonitor.left,
        mi.rcMonitor.bottom - mi.rcMonitor.top,                     // Size and position
        
        NULL,                                                       // Parent window    
        NULL,                                                       // Menu
        hInstance,                                                  // Instance handle
        0                                                           // Additional application data
    );
    
    // Show it
    ShowWindow(hwnd, TRUE);
    UpdateWindow(hwnd);
    
    // Create OpenGL context
    PIXELFORMATDESCRIPTOR pfd =
    {
        sizeof(PIXELFORMATDESCRIPTOR),
        1,
        PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER,    //Flags
        PFD_TYPE_RGBA,        // The kind of framebuffer. RGBA or palette.
        32,                   // Colordepth of the framebuffer.
        0, 0, 0, 0, 0, 0,
        0,
        0,
        0,
        0, 0, 0, 0,
        24,                   // Number of bits for the depthbuffer
        8,                    // Number of bits for the stencilbuffer
        0,                    // Number of Aux buffers in the framebuffer.
        PFD_MAIN_PLANE,
        0,
        0, 0, 0
    };
    
    HDC hdc = GetDC(hwnd);
    
    int  pf = ChoosePixelFormat(hdc, &pfd); 
    SetPixelFormat(hdc, pf, &pfd);
    
    HGLRC glrc = wglCreateContext(hdc);
    wglMakeCurrent (hdc, glrc);
    
    // OpenGL extensions
    glGetProgramiv = (PFNGLGETPROGRAMIVPROC) wglGetProcAddress("glGetProgramiv");
    glGetShaderiv = (PFNGLGETSHADERIVPROC) wglGetProcAddress("glGetShaderiv");
    glGetShaderInfoLog = (PFNGLGETSHADERINFOLOGPROC) wglGetProcAddress("glGetShaderInfoLog");
    glCreateShader = (PFNGLCREATESHADERPROC) wglGetProcAddress("glCreateShader");
    glCreateProgram = (PFNGLCREATEPROGRAMPROC) wglGetProcAddress("glCreateProgram");
    glShaderSource = (PFNGLSHADERSOURCEPROC) wglGetProcAddress("glShaderSource");
    glCompileShader = (PFNGLCOMPILESHADERPROC) wglGetProcAddress("glCompileShader");
    glAttachShader = (PFNGLATTACHSHADERPROC) wglGetProcAddress("glAttachShader");
    glLinkProgram = (PFNGLLINKPROGRAMPROC) wglGetProcAddress("glLinkProgram");
    glUseProgram = (PFNGLUSEPROGRAMPROC) wglGetProcAddress("glUseProgram");
    glGetUniformLocation = (PFNGLGETUNIFORMLOCATIONPROC) wglGetProcAddress("glGetUniformLocation");
    glUniform2f = (PFNGLUNIFORM2FPROC) wglGetProcAddress("glUniform2f");
    glUniform1f = (PFNGLUNIFORM1FPROC) wglGetProcAddress("glUniform1f");
    glGenFramebuffers = (PFNGLGENFRAMEBUFFERSPROC) wglGetProcAddress("glGenFramebuffers");
    glBindFramebuffer = (PFNGLBINDFRAMEBUFFERPROC) wglGetProcAddress("glBindFramebuffer");
    glFramebufferTexture2D = (PFNGLFRAMEBUFFERTEXTURE2DPROC) wglGetProcAddress("glFramebufferTexture2D");
    glNamedRenderbufferStorageEXT = (PFNGLNAMEDRENDERBUFFERSTORAGEEXTPROC) wglGetProcAddress("glNamedRenderbufferStorage");
    glActiveTexture = (PFNGLACTIVETEXTUREPROC) wglGetProcAddress("glActiveTexture");
    glUniform1i = (PFNGLUNIFORM1IPROC) wglGetProcAddress("glUniform1i");
    
    // Browse shaders folder for shaders
    WIN32_FIND_DATA data;
    HANDLE hfile = FindFirstFile(".\\shaders\\*.frag", &data);
    char **filenames = (char **)malloc(sizeof(char*));
    filenames[0] = (char*)malloc(strlen(data.cFileName)+2+strlen(".\\shaders\\"));
    sprintf(filenames[0], ".\\shaders\\%s", data.cFileName);
    printf("Found %s\n", filenames[0]);
    nfiles = 1;
    while(FindNextFile(hfile, &data))
    {
        ++nfiles;
        filenames = (char**)realloc(filenames, nfiles*sizeof(char*));
        filenames[nfiles-1] = (char*)malloc(strlen(data.cFileName)+2+strlen(".\\shaders\\"));
        sprintf(filenames[nfiles-1], ".\\shaders\\%s", data.cFileName);
        printf("Found %s\n", filenames[nfiles-1]);
    } 
    FindClose(hfile);
    printf("Read %d files.\n", nfiles);
    
    // Load shaders
    handles = (int*)malloc(nfiles*sizeof(int));
    programs = (int*)malloc(nfiles*sizeof(int));
    time_locations = (int*)malloc(nfiles*sizeof(int));
    resolution_locations = (int*)malloc(nfiles*sizeof(int));
    scale_locations = (int*)malloc(nfiles*sizeof(int));
    nbeats_locations = (int*)malloc(nfiles*sizeof(int));
    highscale_locations = (int*)malloc(nfiles*sizeof(int));
    fft_texture_locations = (int*)malloc(nfiles*sizeof(int));
    fft_width_locations = (int*)malloc(nfiles*sizeof(int));
    for(int i=0; i<nfiles; ++i)
    {
        printf("Loading Shader %d\n", i);
        
        FILE *f = fopen(filenames[i], "rt");
        if(f == 0)printf("Failed to open file: %s\n", filenames[i]);
        fseek(f, 0, SEEK_END);
        int filesize = ftell(f);
        fseek(f, 0, SEEK_SET);
        char *source = (char*)malloc(filesize+2);
        fread(source, 1, filesize, f);
        fclose(f);
        
        handles[i] = glCreateShader(GL_FRAGMENT_SHADER);
        programs[i] = glCreateProgram();
        glShaderSource(handles[i], 1, (GLchar **)&source, &filesize);
        glCompileShader(handles[i]);
        debug(handles[i]);
        glAttachShader(programs[i], handles[i]);
        glLinkProgram(programs[i]);
        
        glUseProgram(programs[i]);
        time_locations[i] = glGetUniformLocation(programs[i], "iTime");
        resolution_locations[i] = glGetUniformLocation(programs[i], "iResolution");
        scale_locations[i] = glGetUniformLocation(programs[i], "iScale");
        nbeats_locations[i] = glGetUniformLocation(programs[i], "iNBeats");
        highscale_locations[i] = glGetUniformLocation(programs[i], "iHighScale");
        fft_texture_locations[i] = glGetUniformLocation(programs[i], "iFFT");
        fft_width_locations[i] = glGetUniformLocation(programs[i], "iFFTWidth");
    }
    
    glUseProgram(programs[0]);
    glViewport(0, 0, w, h);
    
    // Initialize FFT texture
    fft_texture_size = (int)log2(NFFT)+1;
    printf("fft texture width is: %d\n", fft_texture_size);
    glGenTextures(1, &fft_texture_handle);
    glBindTexture(GL_TEXTURE_2D, fft_texture_handle);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, fft_texture_size, fft_texture_size, 0, GL_RGBA, GL_BYTE, (short*)power_spectrum);
    
    // Set render timer to 60 fps
    UINT_PTR t = SetTimer(hwnd, 1, 1000./60., NULL);
    
    // Get start time for relative time sync
    SYSTEMTIME st_start;
    GetSystemTime(&st_start);
    t_start = (float)st_start.wMinute*60.+(float)st_start.wSecond+(float)st_start.wMilliseconds/1000.;
    
    //FFTW3 Setup
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NFFT);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NFFT);
    p = fftw_plan_dft_1d(NFFT, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    
    // Init sound capture
    WAVEFORMATEX wfx;
    wfx.wFormatTag = WAVE_FORMAT_PCM;
    wfx.nChannels = 1;                    
    wfx.nSamplesPerSec = samplerate;
    wfx.wBitsPerSample = 16;
    wfx.nBlockAlign = wfx.wBitsPerSample * wfx.nChannels / 8;
    wfx.nAvgBytesPerSec = wfx.nBlockAlign * wfx.nSamplesPerSec;
    
    int result = waveInOpen(&wi,            
                WAVE_MAPPER,    
                &wfx,           
                NULL,NULL,      
                CALLBACK_NULL | WAVE_FORMAT_DIRECT  
              );
    printf("WaveInOpen: %d\n", result);
    
    int bsize = buffer_size*wfx.wBitsPerSample*wfx.nChannels/8;
    char * buffers;
    if(double_buffered == 1)
        buffers = (char*)malloc(2*bsize);
    else
        buffers = (char*)malloc(bsize);
    
    for(int i = 0; i < double_buffered+1; ++i)
    {
        printf("Buffer i:\n");
        headers[i].lpData =         buffers+i*bsize;             
        headers[i].dwBufferLength = bsize;
        result = waveInPrepareHeader(wi, &headers[i], sizeof(headers[i]));
        printf("WaveInPrepareHeader: %d\n", result);
        result = waveInAddBuffer(wi, &headers[i], sizeof(headers[i]));
        printf("WaveInAddBuffer: %d\n", result);
    }
    
    result = waveInStart(wi);
    printf("WaveInStart: %d\n", result);
    
    for(int i=0; i<NFFT; ++i)
        values[i] = 0.;
    
    // Main loop
    MSG msg;
    while(GetMessage(&msg, NULL, 0, 0) > 0)
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg); 
    }
    
    return msg.wParam;
}
