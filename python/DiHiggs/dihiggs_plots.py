#!/usr/bin/env python

import ROOT as r
r.gROOT.SetBatch(1)
r.gStyle.SetOptStat(0)

import sys

filedir = "/data/uclhc/uci/user/dantrim/TruthAnalysis/DiHiggs/"
ttbar_file = "%ssuperTruth_ttbar.root"%filedir
x1000_file = "%ssuperTruth_X1000.root"%filedir
x800_file = "%ssuperTruth_X800.root"%filedir
x600_file = "%ssuperTruth_X600.root"%filedir
x400_file = "%ssuperTruth_X400.root"%filedir


class Sample() :
    def __init__(self, filename, name, displayname, color, is_signal) :
        self.filename = filename
        self.name = name
        self.displayname = displayname
        self.color = color
        self.signal = is_signal

        self.tree = None
        self.get_tree(filename) 

        print "Sample %s : %d entries"%(name, self.tree.GetEntries())

    def get_tree(self, filename) :
        tree_name = "superTruth"
        chain = r.TChain(tree_name)
        chain.Add(filename)
        self.tree = chain

def get_variables() :

    variables = {}
    variables["l_pt[0]"] = [5, 0, 300]
    variables["l_pt[1]"] = [5, 0, 300]
    variables["l_eta[0]"] = [0.2, -3, 3]
    variables["l_eta[1]"] = [0.2, -3, 3]
    variables["abs(dphi_ll)"] = [0.1, 0, 3.2]
    variables["pTll"] = [20, 0, 600]
    variables["dRll"] = [0.2, 0, 6]
    variables["abs(cosThetaB)"] = [0.05, 0, 1]
    variables["abs(cosTheta1)"] = [0.05, 0, 1]
    variables["abs(cosTheta2)"] = [0.05, 0, 1]
    variables["abs(dphi_boost_ll)"] = [0.1, 0, 3.2]
    variables["n_jet"] = [1, 0, 12]
    variables["n_ljet"] = [1,0,12]
    variables["n_bjet"] = [1,0,12]
    variables["n_cjet"] = [1,0,12]
    variables["n_hjet"] = [1,0,12]
    variables["n_gjet"] = [1,0,12]
    for i in xrange(3) :
        name = "j_pt[%d]"%i
        variables[name] = [10,0,500]
        name = "j_eta[%d]"%i
        variables[name] = [0.5, -5,5]

        name = "lj_pt[%d]"%i
        variables[name] = [10,0,500]
        name = "lj_eta[%d]"%i
        variables[name] = [0.5, -5,5]
    
        name = "bj_pt[%d]"%i
        variables[name] = [20,0,800]
        name = "bj_eta[%d]"%i
        variables[name] = [0.5, -5,5]

    variables["mass_bb"] = [5,0,350]
    variables["dRbb"] = [0.2, 0, 6]
    variables["dphi_bb"] = [0.1, 0, 3.2]
    variables["abs(dphi_ll_bb)"] = [0.1, 0, 3.2]
    variables["dR_ll_bb"] = [0.2, 0, 6]
    variables["abs(dphi_WW_bb)"] = [0.1, 0, 3.2]
    variables["mass_X"] = [20, 100, 1300]
    variables["mass_X_scaled"] = [20, 100, 1300]
    variables["CMboostMag"] = [0.05, 0, 1]
    variables["abs(dphi_boost_CM_ll_met)"] = [0.1, 0, 3.2]
    variables["met_et"] = [10, 0, 500]
    variables["abs(met_phi)"] = [0.1, 0, 3.2]
    variables["abs(dphi_met_ll)"] = [0.1, 0, 3.2]
    variables["mass_met_ll"] = [10, 0, 400]
    variables["met_pTll"] = [20, 0, 800]
    variables["abs(dphi_boost_ll_met)"] = [0.1, 0, 3.2]

    variables["HT2"] = [20, 0, 1200]
    variables["HT2Ratio"] = [0.05, 0, 1]
    variables["HT2_noMET"] = [20, 0, 1200]
    variables["HT2Ratio_noMET"] = [0.05, 0, 1]
    variables["MTCMS"] = [10, 0, 300]
    variables["MT_HWW"] = [10, 0, 300]
    variables["MT_1"] = [20,0,1200]
    variables["MT_2"] = [10,0,800]
    variables["MT_1_scaled"] = [20,0,1200]
    variables["MT_2_scaled"] = [10,0,800]
    variables["mll"] = [10,0,600]

    return variables

def make_name_ok(var) :
    """
    make names OK for saving
    """
    name = var.replace("abs(","").replace(")","")
    name = name.replace("/","_over_")
    name = name.replace("[","").replace("]","")
    return name

def make_plot(samples, var, bounds) :

    print " > plot %s "%var
    ok_name = make_name_ok(var)

    c = r.TCanvas("c_%s"%ok_name, "", 800, 600)
    c.cd()

    n_bins = bounds[2] - bounds[1]
    n_bins = n_bins / bounds[0]

    #print "NBINS = %d, %d %d"%(n_bins, bounds[1], bounds[2])

    axis = r.TH1F("axis_%s"%ok_name, "", int(n_bins), bounds[1], bounds[2])
    axis.GetXaxis().SetTitle(var)
    axis.GetXaxis().SetTitleFont(42)
    axis.GetXaxis().SetLabelFont(42)

    axis.GetYaxis().SetTitle("Arb. Units")
    axis.GetYaxis().SetTitleFont(42)
    axis.GetYaxis().SetLabelFont(42)

    c.SetGrid(1,1)

    selection = "nLeptons==2 && n_bjet==2"
    #selection = "nLeptons==2 && n_bjet==2 && abs(mass_bb - 125.09) < 20 && abs(dphi_met_ll)<1.5" 
    #selection = "nLeptons==2 && n_bjet==2 && HT2Ratio>0.7 && abs(dphi_met_ll)<1.2 && dRll<1.2 && dRll>0.2"

    leg = r.TLegend(0.7, 0.72, 0.93, 0.93)
    maxy = []

    ttbar_integral = 0.0

    histos = []
    for s in samples :
        h = r.TH1F("h_%s_%s"%(ok_name, s.name), "", int(n_bins), bounds[1], bounds[1])
        h.SetLineColor(s.color)
        h.SetFillColor(0)
        h.SetLineWidth(2)
        if not s.signal :
            #continue
            #h.SetLineStyle(1)
            h.SetFillColorAlpha(2, 0.3)
        h.Sumw2()

        #if s.signal :
        #    selection = "( %s ) * 1e3"%selection
        cut = r.TCut(selection)
        sel = r.TCut("1")

        cmd = "%s>>%s"%(var, h.GetName())
        s.tree.Draw(cmd, cut * sel, "goff")

        #normalize
        err = r.Double(0.0)
        integral = h.IntegralAndError(0,-1, err)
        #if "X" not in s.name :
        #    ttbar_integral = integral
        #else :
        #    integral = ttbar_integral
        if integral != 0 : h.Scale(1/integral)

        leg.AddEntry(h, s.displayname, "l")
        histos.append(h)
        maxy.append(h.GetMaximum())

    
    maxy = 1.25*max(maxy)
    #maxy = 5

    c.cd()
    #c.SetLogy(True)
    axis.SetMaximum(maxy)
    axis.Draw("axis")
    c.Update()

    for h in histos :
        #h.SetMinimum(1e-2)
        h.SetMaximum(maxy)
        h.Draw("hist same")

    leg.Draw()

    text = r.TLatex()
    text.SetNDC()
    text.SetTextFont(42)
    text.DrawLatex(0.15, 0.83, "#bf{#it{ATLAS}} Simulation (Truth)")
    text.DrawLatex(0.15, 0.78, "pp #rightarrow X #rightarrow hh")
    text.DrawLatex(0.15, 0.73, "WWbb selection")
    c.Update()

    outname = ok_name + ".pdf"
    c.SaveAs(outname)

def make_plots(samples) :

    print "make_plots"

    variables = get_variables()

    for var, bounds in variables.iteritems() :
        make_plot(samples, var, bounds)




def main() :
    print "dihiggs plots"

    x1000 = Sample(x1000_file, "X1000", "X (1000 GeV)", 4, True)
    x800 = Sample(x800_file, "X800", "X (800 GeV)", 3, True)
    x600 = Sample(x600_file, "X600", "X (600 GeV)", 6, True)
    x400 = Sample(x400_file, "X400", "X (400 GeV)", 7, True)
    ttbar = Sample(ttbar_file, "ttbar", "t#bar{t}", 1, False)

    samples = [ttbar, x1000, x800, x600, x400]

    make_plots(samples)


#________________________________________________________
if __name__ == "__main__" :
    main()
