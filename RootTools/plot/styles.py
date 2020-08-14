import ROOT

def errorStyle( color, markerStyle = 20, markerSize = 1, width = 1):
    def func( histo ):
        histo.SetLineColor( color )
        histo.SetMarkerSize( markerSize )
        histo.SetMarkerStyle( markerStyle )
        histo.SetMarkerColor( color )
#        histo.SetFillColor( color )
        histo.SetLineWidth( width )
        histo.drawOption = "e0"
        return 
    return func

def invisibleStyle ( ):
    def func( histo ):
        histo.SetMarkerSize(0)
        histo.SetLineWidth(0)
        return
    return func

def lineStyle( color, width = None, dotted=False, dashed=False, errors = False):
    def func( histo ):
        histo.SetLineColor( color )
        histo.SetMarkerSize( 0 )
        histo.SetMarkerStyle( 0 )
        histo.SetMarkerColor( color )
        histo.SetFillColor( 0 )
        if dotted: histo.SetLineStyle( 3 )
        if dashed: histo.SetLineStyle( 7 )
        if width: histo.SetLineWidth( width )
        histo.drawOption = "hist"
        if errors: histo.drawOption+='E'
        return 
    return func

def fillStyle( color, lineColor = ROOT.kBlack, lineWidth = 1, errors = False):
    def func( histo ):
        lc = lineColor if lineColor is not None else color
        histo.SetLineColor( lc )
        histo.SetMarkerSize( 0 )
        histo.SetMarkerStyle( 0 )
        histo.SetMarkerColor( lc )
        histo.SetFillColor( color )
        histo.SetLineWidth( lineWidth )
        histo.drawOption = "hist"
        if errors: histo.drawOption+='E'
        return 
    return func


def styler( **kwargs ):
    """
        kwargs can be things like lineColor or fillColor (i.e. histogram attribites like SetLineColor and SetFillColor etc)
    """
    def func( histo ):
        for k,v in kwargs.items():
            if isinstance(k,str):
                #assume k is something like lineColor and use SetLineColor to style histo
                k_ = f"Set{k[0].upper()+''.join(k[1:])}"
                if hasattr(histo,k_):
                    getattr(histo,k_)(v)
                else:
                    raise NotImplementedError(f"Don't know what to do with style arguments {(k,v)}, amongst {kwargs}....")
                    
            else:
                #assume v is a function
                raise NotImplementedError(f"Don't know what to do with style arguments {(k,v)}, amongst {kwargs}....")
        return
    return func
