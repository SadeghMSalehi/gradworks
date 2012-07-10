#!/usr/bin/python

import io
import Image,ImageOps
from os import curdir, sep
from BaseHTTPServer import BaseHTTPRequestHandler, HTTPServer

class MyHandler(BaseHTTPRequestHandler):
	def image2raw(self, imgfilename):
		im = Image.open(imgfilename)
		if im.size[0] != 256 or im.size[1] != 256:
			im = ImageOps.fit(im, (256,256), Image.ANTIALIAS)
		self.wfile.write(im.tostring())

	def do_GET(self):
		try:
			if self.path.endswith(".jpg") or self.path.endswith(".png") or self.path.endswith(".gif"):
				uri = curdir + sep + self.path
				self.send_response(200)
				self.send_header('Content-type',	'text/plain; charset: x-user-defined')
				self.end_headers()
				self.image2raw(uri)
			else:
				f = open(curdir + sep + self.path) #self.path has /test.html
				self.send_response(200)
				self.send_header('Content-type',	'text/html')
				self.end_headers()
				self.wfile.write(f.read())
				f.close()
		except IOError:
			self.send_error(404,'File Not Found: %s' % self.path)
		return
	
def main():
 try:
	port = 8081
	server = HTTPServer(('', 8081), MyHandler)
	print 'started httpserver... %d' % (port)
	server.serve_forever()
 except KeyboardInterrupt:
	print '^C received, shutting down server'
	server.socket.close()

if __name__ == '__main__':
	main()
