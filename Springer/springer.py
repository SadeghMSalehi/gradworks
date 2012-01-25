#!/usr/bin/python

from pyquery import PyQuery as pq
import sys,os,re

from optparse import OptionParser


__USERAGENT = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/535.7 (KHTML, like Gecko) Chrome/16.0.912.75 Safari/535.7"

def individualDownload(ul, bookIdList, pdflinks, reallyRun):
	cnt = 1
	for k,dataCode in enumerate(bookIdList):
		title = ul("li[data-code=%s]>p.title>a" % (dataCode.upper()))
		href = pdflinks[dataCode]
		if "fulltext" in href:
			fname = (str(cnt) if cnt >= 10 else "0" + str(cnt)) + "_" + title.html()
			fname = fname.replace(" ","_").replace(":","_").replace("&", "and").replace("/","-").replace("?","") + ".pdf"
			cmd = "wget --user-agent=\"%s\" \"%s\" -O \"%s\"" % (__USERAGENT, href, fname)
			print cmd
			if (reallyRun):
				os.system(cmd)
				cnt += 1

def collectionDownload(fout, html, bookIdList, hrefs, run):
	beg = html.find("documentPdfDownloadUrls")
	end = html.find("]", beg + 10)

	if beg > 0:
		if run:
			os.system("rm ???.pdf")
		pdfUrls = html[beg:end+1]
		pdfLinks = re.findall(",*\"(.[a-zA-Z0-9/\-]*\.pdf)\"", pdfUrls)
		print pdfLinks
		cnt = 1
		for url in pdfLinks:
			fname = "%03d.pdf" % cnt
			href = "http://www.springerlink.com/%s" % (url[1:])
			cmd = "wget --user-agent=\"%s\" \"%s\" -O \"%s\"" % (__USERAGENT, href, fname)
			cnt += 1
			print cmd
			if run:
				os.system(cmd)
				pass
		if run:
			joinCMD = "ls ???.pdf | sort | xargs pdfjam --outfile \"%s\"" % fout
			print joinCMD
			os.system(joinCMD)
	else:
		downloaded = {}
		cnt = 1
		print hrefs
		for i,dataCode in enumerate(bookIdList):
			if dataCode not in downloaded:
				downloaded[dataCode] = True
				fname = str(cnt) if cnt > 9 else "0" + str(cnt)
				fname = fname + ".pdf"
				cmd = "wget --user-agent=\"%s\" \"%s\" -O \"%s\"" % (__USERAGENT, hrefs[dataCode], fname)
				print hrefs[dataCode]
				if run:
					os.system(cmd)
					cnt += 1


def main(opts, arg):
	d = pq(filename=opts.html)
	ul = d("ul.manifest")

	pdflinks = {}
	bookIdList = []
	"href extraction"
	for li in ul("li"):
		href = pq(li)("a.pdf-resource-sprite").attr.href
		if href is not None:
			bookId = href.split("/")[2]
			bookHref = "http://www.springerlink.com/%s" % (href[1:])
			if bookId not in pdflinks:
				bookIdList.append(bookId)
			pdflinks[bookId] = bookHref
	"title extract"
	if opts.aggregate:
		collectionDownload(opts.fout, d.html(), bookIdList, pdflinks, opts.run)
	else:
		individualDownload(ul, bookIdList, pdflinks, opts.run)

if (__name__ == "__main__"):
	parser = OptionParser(usage="%prog [options]")
	parser.add_option("-r", "--run", help="run actual download", dest="run", action="store_true", default=False)
	parser.add_option("-a", "--aggregate", help="aggrgate all documents", dest="aggregate", action="store_true", default=False)
	parser.add_option("-i", "--input", help="downloaded html containing pdf links", default="index.html", dest="html")
	parser.add_option("-o", "--output", help="the name of the ouptut file", default="out.pdf", dest="fout")
	(opts, args) = parser.parse_args()
	main(opts, args)
