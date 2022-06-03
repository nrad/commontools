import ROOT as rt

def myGrid( gridOn):
  myStyle.SetPadGridX(gridOn)
  myStyle.SetPadGridY(gridOn)

#fixOverlay: Redraws the axis
def fixOverlay(): gPad.RedrawAxis()

def setMyStyle():
  myStyle =  rt.TStyle("myStyle","Style for P-My")

   #for the canvas:
  myStyle.SetCanvasBorderMode(0)
  myStyle.SetCanvasColor(rt.kWhite)
  myStyle.SetCanvasDefH(600) #Height of canvas
  myStyle.SetCanvasDefW(600) #Width of canvas
  myStyle.SetCanvasDefX(0)   #POsition on screen
  myStyle.SetCanvasDefY(0)


  myStyle.SetPadBorderMode(0)
  #myStyle.SetPadBorderSize(Width_t size = 1)
  myStyle.SetPadColor(rt.kWhite)
  myStyle.SetPadGridX(False)
  myStyle.SetPadGridY(False)
  myStyle.SetGridColor(0)
  myStyle.SetGridStyle(3)
  myStyle.SetGridWidth(1)

#For the frame:
  myStyle.SetFrameBorderMode(0)
  myStyle.SetFrameBorderSize(1)
  myStyle.SetFrameFillColor(0)
  myStyle.SetFrameFillStyle(0)
  myStyle.SetFrameLineColor(1)
  myStyle.SetFrameLineStyle(1)
  myStyle.SetFrameLineWidth(1)
  
#For the histo:
  #myStyle.SetHistFillColor(1)
  #myStyle.SetHistFillStyle(0)
  myStyle.SetHistLineColor(1)
  myStyle.SetHistLineStyle(0)
  myStyle.SetHistLineWidth(1)
  #myStyle.SetLegoInnerR(Float_t rad = 0.5)
  #myStyle.SetNumberContours(Int_t number = 20)

  myStyle.SetEndErrorSize(2)
  #myStyle.SetErrorMarker(20)
  #myStyle.SetErrorX(0.)
  
  myStyle.SetMarkerStyle(20)
  
#For the fit/function:
  myStyle.SetOptFit(1)
  myStyle.SetFitFormat("5.4g")
  myStyle.SetFuncColor(2)
  myStyle.SetFuncStyle(1)
  myStyle.SetFuncWidth(1)

#For the date:
  myStyle.SetOptDate(0)
  # myStyle.SetDateX(Float_t x = 0.01)
  # myStyle.SetDateY(Float_t y = 0.01)

# For the statistics box:
  myStyle.SetOptFile(0)
  myStyle.SetOptStat(0) # To display the mean and RMS:   SetOptStat("mr")
  myStyle.SetStatColor(rt.kWhite)
  myStyle.SetStatFont(42)
  myStyle.SetStatFontSize(0.025)
  myStyle.SetStatTextColor(1)
  myStyle.SetStatFormat("6.4g")
  myStyle.SetStatBorderSize(1)
  myStyle.SetStatH(0.1)
  myStyle.SetStatW(0.15)
  # myStyle.SetStatStyle(Style_t style = 1001)
  # myStyle.SetStatX(Float_t x = 0)
  # myStyle.SetStatY(Float_t y = 0)

# Margins:
  myStyle.SetPadTopMargin(0.05)
  myStyle.SetPadBottomMargin(0.13)
  myStyle.SetPadLeftMargin(0.16)
  myStyle.SetPadRightMargin(0.05)

# For the Global title:

  myStyle.SetOptTitle(0)
  myStyle.SetTitleFont(42)
  myStyle.SetTitleColor(1)
  myStyle.SetTitleTextColor(1)
  myStyle.SetTitleFillColor(10)
  myStyle.SetTitleFontSize(0.05)
  myStyle.SetTitleBorderSize(0)
  myStyle.SetTitleAlign(11)
  myStyle.SetTitleFillColor(0)
  myStyle.SetTitleX(0.5)
  myStyle.SetTitleY(0.95)
  myStyle.SetTitleFont(40)
  myStyle.SetTitleFontSize(0.04)

  # myStyle.SetTitleH(0) # Set the height of the title box
  # myStyle.SetTitleW(0) # Set the width of the title box
  # myStyle.SetTitleX(0) # Set the position of the title box
  # myStyle.SetTitleY(0.985) # Set the position of the title box
  # myStyle.SetTitleStyle(Style_t style = 1001)
  # myStyle.SetTitleBorderSize(2)

# For the axis titles:

  myStyle.SetTitleColor(1, "XYZ")
  myStyle.SetTitleFont(42, "XYZ")
  myStyle.SetTitleSize(0.06, "XYZ")
  # myStyle.SetTitleXSize(Float_t size = 0.02) # Another way to set the size?
  # myStyle.SetTitleYSize(Float_t size = 0.02)
  myStyle.SetTitleXOffset(0.9)
  myStyle.SetTitleYOffset(1.25)
  # myStyle.SetTitleOffset(1.1, "Y") # Another way to set the Offset

# For the axis labels:

  myStyle.SetLabelColor(1, "XYZ")
  myStyle.SetLabelFont(42, "XYZ")
  myStyle.SetLabelOffset(0.007, "XYZ")
  myStyle.SetLabelSize(0.05, "XYZ")

# For the axis:

  myStyle.SetAxisColor(1, "XYZ")
  myStyle.SetStripDecimals(True)
  myStyle.SetTickLength(0.03, "XYZ")
  myStyle.SetNdivisions(510, "XYZ")
  myStyle.SetPadTickX(1)  # To get tick marks on the opposite side of the frame
  myStyle.SetPadTickY(1)

# Change for log plots:
  myStyle.SetOptLogx(0)
  myStyle.SetOptLogy(0)
  myStyle.SetOptLogz(0)

# Postscript options:
  myStyle.SetPaperSize(20.,20.)
  # myStyle.SetLineScalePS(Float_t scale = 3)
  # myStyle.SetLineStyleString(Int_t i, const char* text)
  # myStyle.SetHeaderPS(const char* header)
  # myStyle.SetTitlePS(const char* pstitle)

  # myStyle.SetBarOffset(Float_t baroff = 0.5)
  # myStyle.SetBarWidth(Float_t barwidth = 0.5)
  # myStyle.SetPaintTextFormat(const char* format = "g")
  # myStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0)
  # myStyle.SetTimeOffset(Double_t toffset)
  # myStyle.SetHistMinimumZero(kTRUE)

  myStyle.SetHatchesLineWidth(5)
  myStyle.SetHatchesSpacing(0.05)

  myStyle.cd()


import ROOT

def setExtraStyle():
    ROOT.gStyle.SetPaintTextFormat("0.2f")
    ROOT.gStyle.SetOptTitle(1)
    ROOT.gStyle.SetTitleStyle(0)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetTitleFont(42)
    ROOT.gStyle.SetTitleFontSize(0.035)
    ROOT.gStyle.SetTitleH(0)      # Set the height of the title box
    ROOT.gStyle.SetTitleW(0)      # Set the width of the title box
    ROOT.gStyle.SetTitleX(0.5)   # Set the position of the title box
    ROOT.gStyle.SetTitleY(0.95) # Set the position of the title box
    ROOT.gStyle.SetTitleAlign(11)
    
    ROOT.gStyle.SetPadTickX(0)
    ROOT.gStyle.SetPadTickY(0)
    
    ROOT.gStyle.SetStatX(0.9)
    ROOT.gStyle.SetStatY(0.92)
    ROOT.gStyle.SetStatFontSize(30)
    ROOT.gStyle.SetStatBorderSize(0)
    ROOT.gStyle.SetStatFont(43)
    ROOT.gStyle.SetStatH(0.2)
    ROOT.gStyle.SetStatW(0.4)


def setPalette(name='temp', ncontours=999):
    """Set a color palette from a given RGB list
    stops, red, green and blue should all be lists of the same length
    see set_decent_colors for an example"""
    from array import array
    
    palettes = {
        'gray':{
                'stops' : [0.00, 0.34, 0.61, 0.84, 1.00], 
                'red'   : [1.00, 0.84, 0.61, 0.34, 0.00], 
                'green' : [1.00, 0.84, 0.61, 0.34, 0.00], 
                'blue'  : [1.00, 0.84, 0.61, 0.34, 0.00], 
        },
        'temp':{
                'stops' : [0.00, 0.5,  1.00],
                'red'   : [0.10, 1.00,   0.90], 
                'green' : [0.70, 0.95,   0.30],
                'blue'  : [0.90, 1.00,   0.10], 
        },

    }
    palette=palettes[name]
    
    s = array('d', palette['stops'])
    r = array('d', palette['red'])
    g = array('d', palette['green'])
    b = array('d', palette['blue'])

    npoints = len(s)
    ROOT.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    ROOT.gStyle.SetNumberContours(ncontours)


