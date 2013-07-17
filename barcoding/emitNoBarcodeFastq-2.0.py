import os
import numpy as np
from pbcore.io import FastqWriter, FastqRecord
from pbcore.io.BasH5Reader import *
from pbtools.pbbarcode.BarcodeH5Reader import BarcodeH5Reader

NULL_BARCODE = 'NOBC'


def main( parser ):
    args = parser.parse_args()

    if not os.path.exists( args.outDir ):
        os.makedirs( args.outDir )
        
    emitNoBCFastqs( args.inputFofn, args.barcodeFofn, args.outDir, args.outFile )


def emitNoBCFastqs( inputFofn_filename, barcodeFofn_filename, outDir, outFile ):
    # step through the bas.h5 and barcode.h5 files and emit 
    # reads for each of these. 
    inputFofn = open(inputFofn_filename).read().splitlines()
    barcodeFofn = open(barcodeFofn_filename).read().splitlines()
    outFastq = [] 
        
    for basFile, barcodeFile in zip(inputFofn, barcodeFofn):
        basH5       = BasH5Reader(basFile)
        bcH5        = BarcodeH5Reader(barcodeFile)
        
        msk = -np.in1d( basH5.sequencingZmws, bcH5.bestDS[:,0], assume_unique=True )
        
        for hn in basH5.sequencingZmws[ msk ]:
            zmw = basH5[ hn ]
            if zmw:
                reads = zmw.subreads()
                if any(reads):
                    for read in reads:
                        outFastq.append(FastqRecord(read.readName, 
                                                    read.basecalls(),
                                                    read.QualityValue()))
 
    with FastqWriter("%s/%s.fastq" % (outDir, outFile)) as w:
        for e in outFastq:
                    w.writeRecord(e)
    

if __name__ == '__main__':

    import argparse
    
    parser = argparse.ArgumentParser(prog='emitNoBarcodeFastq', description='write non-barcoded subreads to fastq')
    parser.add_argument('inputFofn', metavar='inputfofn', type=str,
                    help='input.fofn from SMRT Analysis barcode job')
    parser.add_argument('barcodeFofn', metavar='barcodeFofn', type=str,
                    help='barcode.fofn from SMRT Analysis barcode job')
    parser.add_argument('--outDir', metavar='outDir', default=os.getcwd(),
                    help='directory to store output fastq file' )
    parser.add_argument('--outFile', metavar='outFile', default=NULL_BARCODE,
                    help='name of output fastq file' )

    main( parser )
