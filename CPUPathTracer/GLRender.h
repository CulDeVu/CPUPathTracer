#pragma once

#include <windows.h>
#include <thread>
#include <GL/gl.h>

#pragma comment(lib, "opengl32.lib")

//std::thread loopingThread;
GLuint tex;
HWND window;
HDC hDC;
HGLRC hGLRC;
HPALETTE hPalette;

bool running = true;

void init()
{
	/* set viewing projection */
	glMatrixMode(GL_PROJECTION);
	//glFrustum(-0.5F, 0.5F, -0.5F, 0.5F, 1.0F, 3.0F);
	glOrtho(0, 1, 1, 0, 0.01, 100);

	/* position viewer */
	glMatrixMode(GL_MODELVIEW);
	glTranslatef(0.0F, 0.0F, -2.0F);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_TEXTURE_2D);

	glGenTextures(1, &tex);
	glBindTexture(GL_TEXTURE_2D, tex);

	// Black/white checkerboard
	/*float pixels[] = {
		0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f,
		1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f
	};*/
	float* pixels = new float[imageWidth * imageHeight * 3];
	for (int i = 0; i < imageWidth * imageHeight; ++i)
	{
		int x = i % imageHeight;
		int y = i / imageHeight;
		color c = buffer[x][y].normalized();
		pixels[3*i + 0] = c.r;
		pixels[3*i + 1] = c.g;
		pixels[3*i + 2] = c.b;
	}
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 2, 2, 0, GL_RGB, GL_FLOAT, pixels);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
}

void updateTex()
{
	float* pixels = new float[imageWidth * imageHeight * 3];
	for (int i = 0; i < imageWidth * imageHeight; ++i)
	{
		int x = i % imageWidth;
		int y = imageHeight - i / imageWidth - 1;
		color c = buffer[x][y].normalized();
		pixels[3 * i + 0] = c.r;
		pixels[3 * i + 1] = c.g;
		pixels[3 * i + 2] = c.b;
	}
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, imageWidth, imageHeight, 0, GL_RGB, GL_FLOAT, pixels);
	delete[] pixels;
}

void redraw()
{
	/* clear color and depth buffers */
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	updateTex();

	/* draw six faces of a cube */
	glBegin(GL_QUADS);
		glTexCoord2f(0, 0);
		glColor3f(1, 1, 1);
		glVertex2f(0, 0);

		glTexCoord2f(1, 0);
		glColor3f(1, 1, 1);
		glVertex2f(1, 0);

		glTexCoord2f(1, 1);
		glColor3f(1, 1, 1);
		glVertex2f(1, 1);

		glTexCoord2f(0, 1);
		glColor3f(1, 1, 1);
		glVertex2f(0, 1);
	glEnd();

	SwapBuffers(hDC);
}

void resize()
{
	/* set viewport to cover the window */
}

void setupPixelFormat(HDC hDC)
{
	PIXELFORMATDESCRIPTOR pfd = {
		sizeof(PIXELFORMATDESCRIPTOR),  /* size */
		1,                              /* version */
		PFD_SUPPORT_OPENGL |
		PFD_DRAW_TO_WINDOW |
		PFD_DOUBLEBUFFER,               /* support double-buffering */
		PFD_TYPE_RGBA,                  /* color type */
		16,                             /* prefered color depth */
		0, 0, 0, 0, 0, 0,               /* color bits (ignored) */
		0,                              /* no alpha buffer */
		0,                              /* alpha bits (ignored) */
		0,                              /* no accumulation buffer */
		0, 0, 0, 0,                     /* accum bits (ignored) */
		16,                             /* depth buffer */
		0,                              /* no stencil buffer */
		0,                              /* no auxiliary buffers */
		PFD_MAIN_PLANE,                 /* main layer */
		0,                              /* reserved */
		0, 0, 0,                        /* no layer, visible, damage masks */
	};
	int pixelFormat;

	pixelFormat = ChoosePixelFormat(hDC, &pfd);
	if (pixelFormat == 0) {
		MessageBox(WindowFromDC(hDC), "ChoosePixelFormat failed.", "Error",
			MB_ICONERROR | MB_OK);
		exit(1);
	}

	if (SetPixelFormat(hDC, pixelFormat, &pfd) != TRUE) {
		MessageBox(WindowFromDC(hDC), "SetPixelFormat failed.", "Error",
			MB_ICONERROR | MB_OK);
		exit(1);
	}
}

void
setupPalette(HDC hDC)
{
	int pixelFormat = GetPixelFormat(hDC);
	PIXELFORMATDESCRIPTOR pfd;
	LOGPALETTE* pPal;
	int paletteSize;

	DescribePixelFormat(hDC, pixelFormat, sizeof(PIXELFORMATDESCRIPTOR), &pfd);

	if (pfd.dwFlags & PFD_NEED_PALETTE) {
		paletteSize = 1 << pfd.cColorBits;
	}
	else {
		return;
	}

	pPal = (LOGPALETTE*)
		malloc(sizeof(LOGPALETTE) + paletteSize * sizeof(PALETTEENTRY));
	pPal->palVersion = 0x300;
	pPal->palNumEntries = paletteSize;

	/* build a simple RGB color palette */
	{
		int redMask = (1 << pfd.cRedBits) - 1;
		int greenMask = (1 << pfd.cGreenBits) - 1;
		int blueMask = (1 << pfd.cBlueBits) - 1;
		int i;

		for (i = 0; i<paletteSize; ++i) {
			pPal->palPalEntry[i].peRed =
				(((i >> pfd.cRedShift) & redMask) * 255) / redMask;
			pPal->palPalEntry[i].peGreen =
				(((i >> pfd.cGreenShift) & greenMask) * 255) / greenMask;
			pPal->palPalEntry[i].peBlue =
				(((i >> pfd.cBlueShift) & blueMask) * 255) / blueMask;
			pPal->palPalEntry[i].peFlags = 0;
		}
	}

	hPalette = CreatePalette(pPal);
	free(pPal);

	if (hPalette) {
		SelectPalette(hDC, hPalette, FALSE);
		RealizePalette(hDC);
	}
}

LRESULT APIENTRY
WndProc(
HWND hWnd,
UINT message,
WPARAM wParam,
LPARAM lParam)
{
	switch (message) {
	case WM_CREATE:
		/* initialize OpenGL rendering */
		hDC = GetDC(hWnd);
		setupPixelFormat(hDC);
		setupPalette(hDC);
		hGLRC = wglCreateContext(hDC);
		wglMakeCurrent(hDC, hGLRC);
		init();
		return 0;
	case WM_DESTROY:
		/* finish OpenGL rendering */
		if (hGLRC) {
			wglMakeCurrent(NULL, NULL);
			wglDeleteContext(hGLRC);
		}
		if (hPalette) {
			DeleteObject(hPalette);
		}
		ReleaseDC(hWnd, hDC);
		PostQuitMessage(0);
		return 0;
	case WM_SIZE:
		/* track window size changes */
		if (hGLRC) {
			resize();
			return 0;
		}
	case WM_PALETTECHANGED:
		/* realize palette if this is *not* the current window */
		if (hGLRC && hPalette && (HWND)wParam != hWnd) {
			UnrealizeObject(hPalette);
			SelectPalette(hDC, hPalette, FALSE);
			RealizePalette(hDC);
			redraw();
			break;
		}
		break;
	case WM_QUERYNEWPALETTE:
		/* realize palette if this is the current window */
		if (hGLRC && hPalette) {
			UnrealizeObject(hPalette);
			SelectPalette(hDC, hPalette, FALSE);
			RealizePalette(hDC);
			redraw();
			return TRUE;
		}
		break;
	case WM_PAINT:
		{
			PAINTSTRUCT ps;
			BeginPaint(hWnd, &ps);
			if (hGLRC) {
				redraw();
			}
			EndPaint(hWnd, &ps);
			return 0;
		}
		break;

	default:
		break;
	}
	return DefWindowProc(hWnd, message, wParam, lParam);
}

void threadLoop()
{
	MSG msg;
	while (running)
	{
		// Check to see if any messages are waiting in the queue
		while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
		{
			// Translate the message and dispatch it to WindowProc()
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}

		// If the message is WM_QUIT, exit the while loop
		if (msg.message == WM_QUIT)
			break;

		//redraw();
		InvalidateRect(window, NULL, FALSE);

		Sleep(0);
	}
}

int WinMain1()
{
	const char* const myclass = "myclass";
	WNDCLASSEX wndclass = { sizeof(WNDCLASSEX), CS_DBLCLKS, WndProc,
		0, 0, GetModuleHandle(0), LoadIcon(0, IDI_APPLICATION),
		LoadCursor(0, IDC_ARROW), HBRUSH(COLOR_WINDOW + 1),
		0, myclass, LoadIcon(0, IDI_APPLICATION) };
	if (RegisterClassEx(&wndclass))
	{
		window = CreateWindowEx(0, myclass, "title",
			WS_OVERLAPPEDWINDOW, 0, 0,
			imageWidth, imageHeight, 0, 0, GetModuleHandle(0), 0);
		if (window)
		{
			ShowWindow(window, SW_SHOWDEFAULT);
			
			threadLoop();
		}
	}
	return 0;
}

std::thread renderWindow;
void threadedRenderWindow()
{
	renderWindow = std::thread(WinMain1);
}
void exitThreadedRenderWindow()
{
	if (renderWindow.joinable())
		renderWindow.join();

	running = false;
}