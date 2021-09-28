import ROOT
__iColor0=1700
_colors  = []
def Color(iColor=None):
    #from ROOT import ddTColor
    myColor = ROOT.TColor.GetFreeColorIndex()
    my_colors = [
                    #( 56, 105, 177 ), #blue
                    ( 70, 165, 234 ), #blue
                    ( 63, 152, 82 ) , # green
                    ( 211, 93, 96 ) , # red
                    ( 144, 102, 167 ), # purple
                    ( 218, 126, 48 ), # orange


                    (189, 61, 235),
                    (102, 0, 102),
                    (255, 228, 206),
                    (74, 95, 245),
                    (27, 193, 209),
                    (0, 26, 80),
                    (254, 189, 38),
                    (230, 72, 30),
                    (43, 150, 82),
                    (22, 209, 92),
                    (193, 197, 158),
                    (245, 208, 38),
                    (235, 27, 54),

    ]
    
    #indx    = (__iColor - __iColor0)%len(my_colors)
    if iColor != None:
        if iColor>len(my_colors):
            return iColor
        else:
            indx = iColor
        #print('color:', iColor,indx)
    else: 
        indx = Color.index%len(my_colors)
    color = my_colors[indx]
    Color.index+=1
    
    tcolor = ROOT.TColor(myColor, *[x/255. for x in color] )
    _colors.append(tcolor)
    return myColor
Color.index=0 ## hack to keep count how many times function is called

