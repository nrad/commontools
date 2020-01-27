def test2():
    form1 = ROOT.TFormula( 'form1', 'abs(sin(x)/x)' )
    sqroot = ROOT.TF1( 'sqroot', 'x*gaus(0) + [3]*form1', 0, 10 )
    sqroot.SetParameters( 10, 4, 1, 20 )
    c=ROOT.TCanvas()
    sqroot.Draw("l")
    c.Draw()
    return sqroot,c
    