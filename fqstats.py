import sys
import HTSeq
import numpy
from matplotlib import pyplot, ticker
from pandas import DataFrame

fastq = HTSeq.FastqReader(sys.argv[1], 'phred')

def gc_content():
    a = numpy.zeros(101)
    t = numpy.zeros(101)
    g = numpy.zeros(101)
    c = numpy.zeros(101)
    nreads=0
    for i in fastq:
        nreads+=1
        seq = list(str(i))
        pos=0
        for j in seq:
            if j == "A":
                a[pos]+=1
            elif j == "T":
                t[pos]+=1
            elif j == "G":
                g[pos]+=1
            elif j == "C":
                c[pos]+=1
            pos+=1;
    aperc = a / nreads * 100
    tperc = t / nreads * 100
    gperc = g / nreads * 100
    cperc = c / nreads * 100
    fig = pyplot.figure()
    pyplot.plot(aperc)
    pyplot.plot(tperc)
    pyplot.plot(gperc)
    pyplot.plot(cperc)
    pyplot.title("GC content per base in " + str(sys.argv[1]))
    ax=pyplot.subplot(111)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    pyplot.xticks(fontsize="4")
    yminaux = int(min(numpy.concatenate((aperc, tperc, gperc, cperc)))) - 10
    ymin = 0 if yminaux < 0 else yminaux
    ymax = int(max(numpy.concatenate((aperc, tperc, gperc, cperc)))) + 10
    pyplot.yticks(range(ymin, ymax,5), fontsize="6")
    pyplot.grid(color='gray', linestyle='dotted')
    return fig

def qualscore():
    df = DataFrame()
    allvec=[]
    for i in fastq:
        qscore = list(i.qual)
        df=df.append([qscore])
    for j in xrange(0,100):
        allvec.append(list(df[j]))
    fig = pyplot.figure()
    bp = pyplot.boxplot(allvec, 0, '', patch_artist=True)
    pyplot.title("Phred quality score per base in " + str(sys.argv[1]))
    ax=pyplot.subplot(111)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    pyplot.xticks(fontsize="4")
    pyplot.axhspan(0, 20, facecolor='red', alpha=0.1)
    pyplot.axhspan(20, 30, facecolor='yellow', alpha=0.1)
    pyplot.axhspan(30, 40, facecolor='green', alpha=0.1)
    pyplot.yticks(range(0,41), fontsize='6')
    pyplot.xticks(fontsize='6')
    pyplot.grid(color='gray', linestyle='dotted')
    pyplot.setp(bp['boxes'], color='blue', alpha=0.5)
    pyplot.setp(bp['whiskers'], color='black', linestyle='solid', alpha=0.5)
    pyplot.setp(bp['medians'], color='red')
    return fig

def main():
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(str(sys.argv[1]) + ".pdf")
    pdf.savefig(gc_content())
    pdf.savefig(qualscore())
    pdf.close()

if __name__ == "__main__":
    main()
