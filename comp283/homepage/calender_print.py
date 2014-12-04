__author__ = 'joohwi'

import datetime

delta2 = datetime.timedelta(2)
delta7 = datetime.timedelta(7)

monday = datetime.date(2015,1,5)
wedday = monday + delta2
friday = wedday + delta2

k = 0

def print_row(d):
    weekdays = ("Monday", "Wednesday", "Friday")
    print d.strftime("%b %0d"),d.strftime("%a"),


def print_rows(n, m, w, f):
    global k
    nn = n
    color = "#f0f0f0" if n%2 == 0 else "#ffffff"
    print "<tr style='background-color: %s'><td>%s</td><td>"%(color, "%02d"%(k) if n not in (27,28,29,38,6) else "-"),
    print_row(m)
    n += 1
    if n not in (28,29,30,7):
        k += 1
    print "</td><td>desc</td><td><!--links--></td></tr>"
    print "<tr style='background-color: %s'><td>%s</td><td>"%(color, "%02d"%(k) if n not in (27,28,29,38,6) else "-"),
    print_row(w)
    n += 1
    if n not in (28,29,30,40):
        k += 1
    print "</td><td>desc</td><td><!--links--></td></tr>"
    print "<tr style='background-color: %s'><td>%s</td><td>"%(color, "%02d"%(k) if n not in (27,28,29,38,6) else "-"),
    print_row(f)
    n += 1
    print "</td><td>desc</td><td><!--links--></td></tr>"
    if n not in (28,29,30,39):
        k += 1

j = 0
while monday.month < 4 or monday.month == 4 and monday.day < 24:
    print_rows(j, monday, wedday, friday)
    monday = monday + delta7
    wedday = monday + delta2
    friday = wedday + delta2
    j += 3