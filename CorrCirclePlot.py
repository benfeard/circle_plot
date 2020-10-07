#!/usr/bin/env python3
#usage: python CircleDMS.py apoe4.ct e4correlations.dat output.ps

"""
Cross Correlation Circle Plots for Proteins
by Benfeard Williams, II, PhD

This code helps you visualize the cross correlations (both positive and negative)
between residues in your protein. Input data should come from molecular dynamics
simulations or experimental data.

It is recommended to minimize the number of correlations include to produce a cleaner
and easier to understand plot.

Please see generateBaseFile.py to produce an appropriate input file based on FASTA data.

ver. 0.1 - Sat Aug 8, 2020
ver. 0.2 - Wed Oct 7, 2020
"""

import sys
import math

def readCT(ct_file):
    #reads a ct file, !reqires header!
    num,seq,bp = [],[],[]
    linesToRead = int(open(ct_file).readlines()[0].rstrip().split()[0])
    
    for i in open(ct_file).readlines()[1:linesToRead+1]:
        a = i.rstrip().split()
        num.append(int(a[0])),seq.append(str(a[1])),bp.append(int(a[4]))
    return num,seq,bp

def processCorrData(input):
    #Read and store Cross Correlation data where i,j are residue pairs as dictionary
    
    formattedData = {'i':[],'j':[],'corr':[]}
    
    for line in open(input).readlines()[1:]:
        line = line.rstrip().split()
        i = int(line[0])
        j = int(line[1])
        corr = float(line[2])
        formattedData['i'].append(i)
        formattedData['j'].append(j)
        formattedData['corr'].append(corr)
        
    return formattedData

def genCorrelString(corrData):
    """
    correlDat is a dict. i, j, correl are names.
    """
    
    line = ""
    count = 0
    length = len(corrData['i'])
    for num in xrange(length):
        if corrData['j'][num] >= 0.10:
            line += '[{0} {1} 0.00 0.50 0.00 {2}]\n'.format(corrData['i'][num],corrData['j'][num],0)
            count += 1
        if corrData['j'][num] <= -0.04:
            line += '[{0} {1} 0.80 0.10 0.80 {2}]\n'.format(corrData['i'][num],corrData['j'][num],0)
            count += 1
            
    return line, count    

def getScaleFactor(length):
    """
    calculates the scaling factor. Was determined heuristically from a range or RNA sizes in CircleCompair
    """
    
    if length > 76:
        return 74.0875 / float(length) + 0.020367
    else:
        return -1.06*math.log10(length) + 5.44

def makeCircle(num, seq, bp, corrData, offset=1):
    """
    function to make a circleplot from two structures, shannon and shape reactivity and slip status.
    diffColor uses the differential coloring scale based on a slope and intercept of 1 and -1. offset
    delineates nt start position.
    
    returns postscript file lines as a string
    """

    correct = (min(num), max(num))
    missing = []
    extra = []
    bpCorrect = 1
    #correct, missing, extra, bpCorrect = compareRNA(ct_pred, ct_correct)
    
    ct_length = len(seq)
    Coloring = "[0.00 0.00 0.00]\n" * ct_length
    
    pairingTable = "[1 " + str(ct_length) + " 0.80 0.80 0.80 0]\n" #numbers are the color
    numPairs = 0 #might need to be set equal to 1
    
    pairingTableCorrel, numPairsCorrel = genCorrelString(corrData)
    
    # add the new correlations to the pairing array
    pairingTable += pairingTableCorrel
    numPairs += numPairsCorrel
    #print pairingTable
    scaleFactor = getScaleFactor(ct_length)
    
    # heuristically based on circlecompare
    cirRadius = 3 * ct_length + 14
    
    sens = bpCorrect/float(len(correct)+len(missing))
    ppv = bpCorrect/float(len(correct)+len(extra))


    Circle = """%!

% Set font size and type.
/fontSize 24 def
/quarterFont fontSize 4 div def
/halfFont fontSize 2 div def
/Courier findfont fontSize scalefont setfont
% Set variables handling number, placement of nucleotides.
/numBases {0} def
/currentBase 0 def
/basePoints numBases array def

% Set variables handling scaling, translation of circular backbone.
/scaleFactor {1} def
/scaleFactorX scaleFactor def
/scaleFactorY scaleFactor def
/translateFactorX 0 def
/translateFactorY 0 def

% Set variables handling properties of circular backbone.
/radius {2} def
/center 306 scaleFactor div def
/angle 360 numBases 2 add div -1 mul def
/labelSpace 40 def

% Create the array of nucleotides.
/bases [
""".format(ct_length,scaleFactor,cirRadius)

    Circle += ' '.join(map('({0})'.format, seq)) # lists sequence as '(A) (G)...'
    Circle +="""
] def

% Write the array of pairings.
/pairings [ """
    Circle += pairingTable
    Circle += """
] def

% Set variables handling number, placement of pairings.
/numPairings {0} def
/numPseudoknotted 0 def
/currentPairing 0 def

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write out the combined circular structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the font size and type for the descriptor.
/Courier-Bold findfont fontSize 2 div scalefont setfont

/numDescriptors 7 def
/descriptors [
] def

% Write the contents of the descriptor array.
descriptors {{ aload pop moveto setrgbcolor show }} forall

% Write the array of annotation colors.
/annotations [ """.format(numPairs, sys.argv[1], sys.argv[1])
    Circle += Coloring
    Circle += """] def

% Create the legend array.
/legend [
    [(CorrCoef >= 0.70) 0.0 0.5 0.0]
    [(CorrCoef <= -0.40) 0.8 0.1 0.8]
] def

% Show the legend.
/Courier-Bold findfont 8 scalefont setfont
40 595 moveto

legend {{ gsave aload pop setrgbcolor show grestore 0 -10 rmoveto }} forall

% Reset the font to its original size, type, and color.
/Courier findfont fontSize scalefont setfont
0 setgray

% Write the translation and scaling of the main image.
gsave
scaleFactorX scaleFactorY scale
translateFactorX translateFactorY translate

% Use repeat loop to write nucleotides on the circular path.
0 1 numBases {{
    clear

    % Move to appropriate point, angle to write next nucleotide.
    center center moveto
    currentBase angle mul rotate
    quarterFont radius rmoveto

    % Determine x and y coordinates of the nucleotide location.
    currentBase angle mul -1 mul rotate
    currentpoint
    /y exch cvi def
    /x exch cvi def
    currentBase angle mul rotate
    quarterFont -1 mul 0 rmoveto

    % Save the nucleotide location and show the nucleotide.
    /point [ x y ] def
    basePoints currentBase point put
    
    %define offset for nums
    /offset {2} def

    % Set variables for conditions where number labels are found.
    /numCond1 currentBase offset add 10 mod 0 eq def
    /numCond2 currentBase offset add 1 eq def
    /numCond3 currentBase offset add numBases eq def

    /annotation annotations currentBase get def
    annotation 0 get annotation 1 get annotation 2 get setrgbcolor

    % If a condition is met for a number label, write a label.
    numCond1 numCond2 or numCond3 or {{
        /Courier-Bold findfont fontSize scalefont setfont
        bases currentBase get show
        halfFont -1 mul labelSpace rmoveto
        /numString 7 string def
        currentBase offset add numString cvs show
        0 labelSpace -1 mul rmoveto
        /Courier findfont fontSize scalefont setfont
    }}
    {{ bases currentBase get show }} ifelse
    0 setgray

    % Return to circle center, rotate to 0, increment base.
    0 radius -1 mul rmoveto
    currentBase angle mul -1 mul rotate
    /currentBase currentBase 1 add def
}} repeat

% Write the pairing lines inside the backbone.
0 setgray
0 1 numPairings {{
    /pair pairings currentPairing get def

    % Determine coordinates of first nucleotide in pair.
    /base1 pair 0 get 1 sub def
    /point1 basePoints base1 get def
    /x1 point1 0 get def
    /y1 point1 1 get def

    % Determine coordinates of second nucleotide in pair.
    /base2 pair 1 get 1 sub def
    /point2 basePoints base2 get def
    /x2 point2 0 get def
    /y2 point2 1 get def

    % If the pair should be colored, set its appropriate color.
    pair 2 get pair 3 get pair 4 get setrgbcolor
    
    % if the pair is slipped set variable so
    /isSlipped pair 5 get def

    % Draw current pair, then increment current pairing.
    /between base2 base1 sub def
    between numBases 2 div gt {{ /between numBases between sub def }} if

    /midX x1 x2 add 2 div def
    /midY y1 y2 add 2 div def

    /gamma 0.9 def
    /centerThresh radius 2 mul 8 div 5 mul def

    /ends x2 x1 sub x2 x1 sub mul y2 y1 sub y2 y1 sub mul add sqrt def
    /distance between 2 mul numBases div radius mul gamma mul def
    /lineAngle center midY sub center midX sub atan def
    /distX lineAngle cos distance mul 2 mul def
    /distY lineAngle sin distance mul 2 mul def

    /controlX midX distX add def
    /controlY midY distY add def

    ends centerThresh ge {{
        /controlX center def
        /controlY center def
    }} if

    isSlipped 0.5 gt {{
    x1 y1 moveto x1 y1 controlX controlY x2 y2 curveto [15 5] 0 setdash stroke
    }}{{
    x1 y1 moveto x1 y1 controlX controlY x2 y2 curveto [100000 1] 0 setdash stroke
    }} ifelse
    0 setgray
    /currentPairing currentPairing 1 add def
}} repeat


 % Set variables handling number, placement of nucleotides.
 /numBasesf {1}""".format(numPairs,ct_length,offset)
    Circle += """ def
 /currentBasef 0 def
 /basePointsf numBasesf array def
 
 
 % Set variables handling properties of circular backbone.
 /radiusf {0} def
 /centerf 306 scaleFactor div def
 /anglef 359.75 numBases 2 add div -1 mul def
 /labelSpace 20 def
 
 /annotationsf [ """.format(int(cirRadius+10)) + Coloring + """
 ] def
 
 
 /basesf [ """
 
    Circle += ' '.join(map('()'.format, seq))
    Circle += """

 ] def
 
 
 
 
 
 /fontSizef 40 def
 /quarterFontf fontSizef 4 div def
 
 
 % Use repeat loop to write nucleotides on the circular path.
 0 1 numBasesf {{
     clear
 
     % Move to appropriate point, angle to write next nucleotide.
     center centerf moveto
     currentBasef anglef mul rotate
     quarterFontf radiusf rmoveto
 
     % Determine x and y coordinates of the nucleotide location.
     currentBasef anglef mul -1 mul rotate
     currentpoint
     /xf exch cvi def
     /yf exch cvi def
     currentBasef angle mul rotate
     quarterFontf -1 mul 0 rmoveto
 
     % Save the nucleotide location and show the nucleotide.
     /point [ xf yf ] def
     basePointsf currentBasef point put
     /offset {5} def
 
     % Set variables for conditions where number labels are found.
     /numCond1 currentBasef offset add 10 mod 0 eq def
     /numCond2 currentBasef offset add 1 eq def
     /numCond3 currentBasef offset add numBases eq def
 
     /annotationf annotationsf currentBasef get def
     annotationf 0 get annotationf 1 get annotationf 2 get setrgbcolor
 
     % If a condition is met for a number label, write a label.
     numCond1 numCond2 or numCond3 or {{
         /Courier-Bold findfont fontSizef scalefont setfont
         basesf currentBasef get show
         /Courier-Bold findfont fontSize scalefont setfont
         halfFont -1 mul labelSpace rmoveto
         /numString 7 string def
        % currentBasef offset add numString cvs show
         0 labelSpace -1 mul rmoveto
         /Courier findfont fontSizef scalefont setfont
     }}
     {{ basesf currentBasef get show }} ifelse
     0 setgray
 
     % Return to circle center, rotate to 0, increment base.
     0 radiusf -1 mul rmoveto
     currentBasef angle mul -1 mul rotate
     /currentBasef currentBasef 1 add def
 }} repeat
 
  % Write the array of annotation colors.
 

grestore

/Courier-Bold findfont fontSize 2 div scalefont setfont

/x1 570 (Accepted:) stringwidth pop sub def
/x2 570 (Pairs: {2}) stringwidth pop sub def
/x3 570 (Pseudoknotted Pairs: 0) stringwidth pop sub def
/x4 570 (Sensitivity: {0} / {2} = {3:.2%}) stringwidth pop sub def
/x5 570 (PPV: {0} / {1} = {4:.2%}) stringwidth pop sub def

/statistics [

] def

statistics {{ aload pop moveto show }} forall
showpage
""".format(len(correct), len(correct)+len(extra), len(correct)+len(missing), sens, ppv, offset)

    return Circle


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print 'Usage: CorrCirclePlot.py structureFile corrData outputPlot.ps'
        sys.exit()

    num, seq, bp = readCT(sys.argv[1])

    corrData = processCorrData(sys.argv[2]) #done and updated
    
    circleLines = makeCircle(num, seq, bp, corrData, offset=1) #look into the two ct_pred
    
    outputPlot = open(sys.argv[3],'w')
    outputPlot.write(circleLines)
    outputPlot.close()
