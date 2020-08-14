import RootTools.core.standard as RootTools
import param
import os


class PlotVariable( param.Parameterized ):
    var    = param.String()
    bins   = param.List()
    name   = param.String()
    xTitle = param.String()
    overflow = param.Selector( default="both", objects=["upper","lower","both", None] )
    cut    = param.String()
    logY   = param.Boolean(default=True)
    
    
    def __init__(self,**kwargs):
        if not 'name' in kwargs:
            kwargs['name']=str(kwargs['var']) ## should probably remove special characters
        super().__init__(**kwargs)
        if not self.xTitle:
            self.xTitle = self.name






def makePlotFromPlotVar(stack, plot_var, selection_string ):
    plot = RootTools.Plot( name = plot_var.name, attribute = plot_var.var, 
                           texX = plot_var.xTitle, addOverFlowBin = plot_var.overflow, 
                           binning = plot_var.bins, selectionString = selection_string, stack=stack,
                         )
    #plot.params = plot_var
    for attr,value in plot_var.get_param_values():
        if not hasattr(plot,attr):
            setattr(plot,attr,value)

    return plot








###
###  RootTools plotting
###

def makeAllPlots(stack, selection, plot_names, plot_variables=None):
    
    pvars = [ plot_variables[p] for p in plot_names] if plot_variables else plot_names
    plots = list( map( lambda pvar: makePlotFromPlotVar( stack, pvar, selection) , pvars ))
    return plots



def drawPlot( plot, legend=([0.2,0.6,0.9,0.9],3), save_dir = "./figs/", sorting = True, logY = True , **kwargs):
    """
        wrapper for RootTools plotting...
    """
    c = RootTools.plotting.draw(  plot, 
                              legend = legend,
                              #plot_directory = "./figs/", 
                              plot_directory=save_dir,
                              sorting = sorting,
                              copyIndexPHP = True,
                              logY = logY,
                              **kwargs,
                            #canvasModifications=[ lambda c: c.Draw()]
                            
                          )
    return c


def plotAndDraw(stack, 
                region, 
                plot_names,
                legend=([0.2,0.6,0.9,0.9],3),
                save_dir_base = "./figs/",
                plot_variables = None,
                **kwargs,
               ):
    """
        wrapper for RootTools plotting
    """
    if not plot_variables:
        if not all( [isinstance(p,PlotVariable) for p in plot_names]):
            raise Exception("If plot_variables is not given, the plot_names should be PlotVariables")

    plots = makeAllPlots(stack, region.formatted, plot_names , plot_variables)
    save_dir = os.path.join( save_dir_base, region.name )
    canvs = []
    for p in plots:
        RootTools.plotting.fill_with_draw([p])
        c = drawPlot(p, legend, save_dir, logY=getattr(p,'logY',True) , **kwargs)
        #c[0].Draw()
        #print(c)
        canvs.append(c)
    return plots, canvs
