# SurfaceGenie_0.1/ui.R
library(shiny)
library(plotly)

scores <- c("GenieScore", "IsoGenie", "OmniGenie", "IsoOmniGenie")
images <- c("gs.png", "isg.png", "og.png", "iog.png")

shinyUI(navbarPage("  SurfaceGenie  ", theme = "bootstrap.css",
  
  ##########    Home    ##########
  
  tabPanel(
    "Home",
    div(
     h4("Welcome to ", span(class ="text-success", "SurfaceGenie"),"!"),
     p(tags$i("Integrating predictive and empirical data for rational marker prioritization")),
      tags$img(src="website_homepage.png", width="330px", align="right"),
     p("SurfaceGenie is a web app for analyzing omic datasets (e.g. proteomic, transcriptomic) to prioritize candidate cell-type specific markers of interest for immunophenotyping, immunotherapy, drug targeting, and other applications. It works by calculating the likelihood a molecule is informative for distinguishing among sample groups (e.g. cell types, experimental conditions). While a major benefit of SurfaceGenie is the ability to prioritize proteins that are localized to the cell surface, it is also possible to analyze data without this parameter to find proteins of interest that reside in other subcellular localizations. See the descriptions for each of the four permutations of the scoring algorithm."),
     p("SurfaceGenie works well with quantitative proteomic and transcriptomic datasets, but others are also possible (see below). All calculations performed within SurfaceGenie are context-dependent, meaning that the tools will consider all data within a single dataset input (which may contain multiple experiments and/or cell types). If a user performs a comparison and subsequently determines additional data should be considered, a new file containing all data for the new comparison is required."),
     br(),
     tags$b("Scoring Permutations"),
     p( tags$i("GenieScore: "), "Use to prioritize ", tags$b("surface proteins"), "that have ", tags$b("disparate"), " levels of abundance/expression."),
     p( tags$i("IsoGenieScore: "), "Use to prioritize ", tags$b("surface proteins"), "that have ", tags$b("similar, high"), " levels of abundance/expression."),
     p( tags$i("OmniGenieScore: "), "Use to prioritize ", tags$b("any molecules"), " (genes/proteins) that have ", tags$b("disparate"), " levels of abundance/expression."),
     p( tags$i("IsoOmniGenieScore: "), "Use to prioritize ", tags$b("any molecules"), " (genes/proteins) that have ", tags$b("similar, high"), " levels of abundance/expression.")
    ),
    br(),
    div(
      tags$b("Overview of Inputs and Outputs"),
      p(tags$i("Input:  ")),
     p(style="margin-left:1.5em", "SurfaceGenie accepts a .csv file containing a list of proteins (UniProt Accession) and a surrogate value representative of abundance (e.g. number of peptide spectrum matches, peak area, FKPM, RKPM) identified within a set of samples. There is no limit to the number of samples that can be analyzed in a single file. SurfaceGenie has SPC datasets for human, mouse, and rat."),
     p(tags$i("Data Processing:  ")),
     p(style="margin-left:1.5em", "SurfaceGenie calculates the dot product of three independent scores:"),
     tags$ol(
       tags$li(tags$u("Surface Protein Consensus (SPC) score"), br(), "A predictive measure of the likelihood 
               that a particular protein can be present at the cell surface. This value is a sum 
               of the number of predictive datasets for which a protein has been predicted to be 
               localized to the cell surface. Scores range 0-4. For more details on the predictive 
               datasets used, see the References tab."),
       tags$li(tags$u("Distribution Score"), br(), "A measure of how evenly or unevenly distributed a 
               protein is among multiple samples within a comparison dataset. It is based on the 
               Gini coefficient for calculating statistical dispersion of values. Scores range 0 - 1/(1-N)."),
       tags$li(tags$u("Signal Strength"), br(), "An approximate measure of protein abundance for cell types 
               in which a protein is observed. Proteins at the lower limit of detection are of lower 
               priority than those with more observations, because it is expected that those of higher 
               abundance will practically serve as more accessible markers for downstream technologies. 
               Scores typically range 0 ~ 4 .")
     ),
     p(tags$i("Output:  ")),
     tags$ul(
       tags$li(tags$u("CSV Download"), br(), "Columns of selected data types (e.g. SPC score, CD molecule annotation, etc) are appended to each entry in the original input file"),
       tags$li(tags$u("Plots"), br(), "Scores from each of the 4 permutations are plotted in order of priority for all proteins within a dataset"),
       tags$li(tags$u("SPC Histogram"), br(), "Displays the distribution of SPC scores")
     ), 
     br(),
     h5(class="text-info", "SPC Score Lookup"),
     p("This feature enables users to obtain Surface Protein Consensus (SPC) score for proteins 
       of interest without analyzing data through SurfaceGenie. Users may perform a batch retrieval 
       by uploading a .csv file containing UniProt Accession numbers or may search individual 
       UniProt accession numbers."),
     br(),
     h5(class="text-info", "Other Applications"),
     p("Although the calculation of SPC score depends on the use of Uniprot Accession IDs (for human, mouse, or rat), the other terms used here are agnostic to the type and distribution of data. Therefore, the OmniGenieScore and IsoOmniGenieScore can be used for any type of quantitative data for which there is a desire to find measurements that are either unique or similar between all samples. This could include metabolomic data, glycomics data, or strain counts for microbiome studies.")
     
    )
  ),

  
  ##########  Instructions ##########

  tabPanel(
    "Instructions",
    div(
      h4("Data Upload Instructions"),
      h5(class="text-info", "Data Format"),
      p("A .csv file containing a list of proteins (UniProt Accession) and a surrogate value 
        representative of abundance (e.g. number of peptide spectrum matches, peak area) 
        identified within a set of samples. "),
      p("The first column of your data file ", tags$b("must be labeled 'Accession'"), " with no extra characters 
        (e.g. not 'Accession #'). This column should contain the UniProt accession numbers of the 
        proteins in your samples. You may include isoforms. To convert from a different protein ID 
        type to UniProt, bulk conversion is available ", a(href="https://www.uniprot.org/uploadlists/", "here"), 
        ". Under 'Select options', select your ID type in the 'From' field and then 'UniProt KB'
        in the 'To' field."
      ),
      p("Importantly, data files ", tags$b("must be in .csv format"), ". If you are working in Excel, click 
        'File --> Save As' and select csv in the drop-down menu to convert from .xlsx to .csv."),
      h5(class="text-info", "Example Data"),
      tableOutput("example_data"),
      br()
      ),
    div(
      h4("Data Proccessing Options"),
      p( tags$i("GenieScore: "), "Use to prioritize ", tags$b("surface proteins"), "that have ", tags$b("disparate"), " levels of abundance/expression."),
      p( tags$i("IsoGenieScore: "), "Use to prioritize ", tags$b("surface proteins"), "that have ", tags$b("similar, high"), " levels of abundance/expression."),
      p( tags$i("OmniGenieScore: "), "Use to prioritize ", tags$b("any molecules"), " (genes/proteins) that have ", tags$b("disparate"), " levels of abundance/expression."),
      p( tags$i("IsoOmniGenieScore: "), "Use to prioritize ", tags$b("any molecules"), " (genes/proteins) that have ", tags$b("similar, high"), " levels of abundance/expression."),
      
      h5(class="text-info", "Sample Grouping"),
      p("Ideally, similar samples such as technical replicates or biological replicates will have
        values averaged or summed together into a single column. However,
        SurfaceGenie will carry out this step for you if you select 'Group samples'. If this box is
        checked, you will need to provide the grouping method as well as the column numbers for each group. 
        For example, If columns 2, 3, and 5 of your dataset should be grouped together and columns 4 & 6 
        comprise another group, you should indicate the presence of 2 groups using the slider and 
        then enter the corresponding column numbers below separated by commas: Group 1: '2, 3, 5', 
        Group 2: '4, 6'. Remember that column 1 will contain accession numbers and cannot be grouped
        with other columns."),
      br()
      ),
    div(
      h4("Data Export Options"),
      h5(class="text-info", "CSV File"),
      p("You may select data to export as columns appended to the right of your original data. 
                        The following variables are available for export:"),
      p("GenieScore Components:"),
      tags$ul(
        tags$li(tags$u("Surface Protein Concensus score (SPC)") ,": A predictive measure of the likelihood that a 
                                particular protein can be present at the cell surface."),
        tags$li(tags$u("Distribution Score (Gini)"), ": A measure of the distribution of the protein amongst samples. 
                                A higher value corresponds to a more localized distribution.",
                a(href="https://en.wikipedia.org/wiki/Gini_coefficient", "Wikipedia - Gini coefficient")),
        tags$li(tags$u("Signal strength (SS)"), ": A weighted value of the maximum value reported among samples for the protein"),
        tags$li(tags$u("Genie Score (GS)"), ": SurfaceGenie's measure for the value of a protein as a potential 
                                marker of interest.")
      ),
      p("Annotations/Linkouts:"),
      tags$ul(
        tags$li(tags$u("CD molecules (CD)"), ": Cluster of differentiation (CD) molecules are annotated with CD nomenclature. CD molecules have validated antibodies against them and therefore are attractive candidate markers for immunodetection -based applications."),
        tags$li(tags$u("HLA molecules (HLA)"), "Human leukocyte antigen (HLA) molecules are surface proteins that have high sequence similarity. As such, it is often challenging to be certain of the specific gene product based solely on peptide-level evidence particularly for Cell Surface Capture experiments. As a result, it may be useful to exclude these from consideration when attempting to identify cell surface makers for a specific cell type."),
        tags$li(tags$u("Number of CSPA experiments (CSPA)"), ": The number of cell types in which this protein was observed in the", tags$a(href="http://wlab.ethz.ch/cspa/", "Cell Surface Protein Atlas."),  "This information can provide context for how specific a protein might be among cell types."),
        tags$li(tags$u("UniProt Linkout"), ": Link to the UniProt entry for input proteins providing effortless access to additional information about candidate markers.")
      ),
      h5(class="text-info", "Plots"),
      p("Several visualizations are made available by SurfaceGenie:"),
      tags$ul(
        tags$li(tags$u("SPC Histogram"), ": Shows the distribution of SPC scores."),
        tags$li(tags$u("Plots"), ": Scores from each of the 4 permutations are plotted in order of priority for all proteins within a dataset.")
      )
    )
    ),
  
  ##########  SurfaceGenie ##########
  
  tabPanel(
    "SurfaceGenie",
    sidebarPanel(
      h5(class="text-info", "Data Input"),
      fileInput("file1", "Choose Input File", multiple =FALSE, 
                accept=c(".csv", ".tsv", ".txt", ".tab", ".xls", ".xlsx"), 
                buttonLabel = "Browse...", placeholder = "No file selected"),
      h5(class="text-info", "Scoring Options"),
      #### regular GS
      #### GS, no SPC = iGenie
      #### GS reversed = eineG
      #### GS reversed, no SPC = eineGi
      
      
      checkboxGroupInput(
        "scoring_opts", label=NULL,
        choiceNames = mapply(scores, images, FUN=function(score, imgloc) {
          tagList(
            score,
            tags$img(src=imgloc, width=75)
          )
        }, SIMPLIFY = FALSE, USE.NAMES = FALSE),
        choiceValues = list(
          "GS", "eineG", "iGenie", "eineGi"),
          selected = list("GS")
      ),

      h1(),
      h5(class="text-info", "Species"),
      radioButtons(
        "species", NULL,
        choices = list(
          "human",
          "rat",
          "mouse"),
        selected = list("human")
        #        choiceNames = NULL,
        #        choiceValues = NULL
      ),      

      h1(),
      h5(class="text-info", "Processing Option"),
      checkboxGroupInput(
        "processing_opts", NULL,
        choiceNames = list(
          "Group samples"),
        choiceValues = list(
          "grouping")
      ),
      conditionalPanel(
        condition = "input.processing_opts.indexOf('smarker') > -1",
        h5(class="text-info", "Markers for Specific Sample"),
        textInput(
          "markersample", "Enter sample name:", placeholder="i.e. 'd00' or 'Group 1'"
        )
      ),
      conditionalPanel(
        condition = "input.processing_opts.indexOf('grouping') > -1",
        h5(class="text-info", "Sample Grouping"),
        selectInput("groupmethod", "Grouping method",
                    choices = c(
                      "Mean" = "ave",
                      "Median" = "med"),
                    selected = "ave"
        ),
        p("*Please see Sample Grouping section on the Home page for instructions 
          on how to enter grouping information."),
        sliderInput("numgroups", "Number of groups",
                    min=2, max=5, value=2, step=1, ticks=FALSE),
        textInput("group1", "Group 1", placeholder="Columns in Group 1"),
        textInput("group2", "Group 2", placeholder="Columns in Group 2"),
        conditionalPanel(
          condition = "input.numgroups >= 3",
          textInput("group3", "Group 3", placeholder="Columns in Group 3")
        ),
        conditionalPanel(
          condition = "input.numgroups >= 4",
          textInput("group4", "Group 4", placeholder="Columns in Group 4")
        ),
        conditionalPanel(
          condition = "input.numgroups >= 5",
          textInput("group5", "Group 5", placeholder="Columns in Group 5")
        ),
        conditionalPanel(
          condition = "input.numgroups >= 6",
          textInput("group6", "Group 6", placeholder="Columns in Group 6")
        ),
        conditionalPanel(
          condition = "input.numgroups >= 7",
          textInput("group7", "Group 7", placeholder="Columns in Group 7")
        ),
        conditionalPanel(
          condition = "input.numgroups >= 8",
          textInput("group8", "Group 8", placeholder="Columns in Group 8")
        ),
        conditionalPanel(
          condition = "input.numgroups >= 9",
          textInput("group9", "Group 9", placeholder="Columns in Group 9")
        ),
        conditionalPanel(
          condition = "input.numgroups >= 10",
          textInput("group10", "Group 10", placeholder="Columns in Group 10")
        )
      ),

      h1(),
      h5(class="text-info", "Export Options (CSV Download Tab)"),
      checkboxGroupInput(
        'export_options1', "SurfaceGenie Components:",
        choiceNames = list(
          "SPC score (SPC)",
          "Gini coefficient (Gini)",
          "Signal strength (SS)"
        ),
        choiceValues = list(
          "SPC", "Gini", "SS"),
        selected = list("GS")
      ),
      checkboxGroupInput(
        'export_options2', "Annotations / Link outs:",
        choiceNames = list(
             "HLA molecules",
             "CD molecules",
             "Gene Name",
             "Number of CSPA experiments",
             "UniProt Linkout", 
             "Transmembrane",
             "Subcellular Location"),
        choiceValues = list(
          "HLA", "CD", "geneName", "CSPA..e", "UniProt Linkout", "Transmembrane", "CC")
  )
),
    
    mainPanel(
      span(textOutput("txtWarning"), style="color:red"),
      tabsetPanel(type = "tabs",
                  tabPanel(
                    "Data Input",
                    tableOutput("data_input"),
                    em(textOutput("input_size"))
                  ),
                  tabPanel(
                    "CSV Download", 
                    tableOutput("data_output"),
                    conditionalPanel(
                      condition = "output.data_output != undefined",
                      em(textOutput("output_size"))
                    ),
                    br(),
                    uiOutput("csv_dlbutton", class="download_this" ), uiOutput("tsv_dlbutton", class="download_this"), uiOutput("xlsx_dlbutton", class="download_this")
                  ),
                  tabPanel(
                    "Plots",
                    plotOutput("SG_SPC_hist"),
                    div(class="bnav",
                    uiOutput("SG_SPC_hist_PNGdlbutton", class="download_this"),
                    uiOutput("SG_SPC_hist_SVGdlbutton", class="download_this")),
                    p(),
                    br(),
                    br(),
                    br(),
                    conditionalPanel(
                      condition = "input.scoring_opts.indexOf('GS')>-1",
                      plotlyOutput("SG_dist"),
                      div(class="bnav",
                      uiOutput("SG_dist_PNGdlbutton", class="download_this"),
                      uiOutput("SG_dist_SVGdlbutton", class="download_this")),
                      p(),
                      br(),
                      br(),
                      br()
                    ),
                    #eineG
                    conditionalPanel(
                      condition = "input.scoring_opts.indexOf('eineG')>-1",
                      plotlyOutput("eineG_dist"),
                      div(class="bnav",
                      uiOutput("eineG_dist_PNGdlbutton", class="download_this"),
                      uiOutput("eineG_dist_SVGdlbutton", class="download_this")),
                      p(),
                      br(),
                      br(),
                      br()
                    ),
                    #iGenie
                    conditionalPanel(
                      condition = "input.scoring_opts.indexOf('iGenie')>-1",
                      plotlyOutput("iGenie_dist"),
                      div(class="bnav",
                      uiOutput("iGenie_dist_PNGdlbutton", class="download_this"),
                      uiOutput("iGenie_dist_SVGdlbutton", class="download_this")),
                      p(),
                      br(),
                      br(),
                      br()
                    ),
                    #eineGi
                    conditionalPanel(
                      condition = "input.scoring_opts.indexOf('eineGi')>-1",
                      plotlyOutput("eineGi_dist"),
                      div(class="bnav",
                      uiOutput("eineGi_dist_PNGdlbutton", class="download_this"),
                      uiOutput("eineGi_dist_SVGdlbutton", class="download_this")),
                      br()
                    )
                  )
      )
    )
  ),
  
  ##########  Surface Protein Concensus (SPC) Score  Lookup ##########
  
  tabPanel(
    "SPC Score  Lookup",
    sidebarPanel(
      h5(class="text-info", "Quick Lookup"),
      textAreaInput("quicklookup", "Uniprot accession number:", 
                placeholder="Enter accession numbers, each on a new line. For example:
                                                              
                                        A0AVT1-1 
                                        A0FGR8-6 
                                        A1L0T0 
                                        A1X283",
                rows=10),
      br(),
      h5(class="text-info", "Bulk Lookup"),
      fileInput("file2", "Choose Input File", multiple =FALSE, 
                accept=c(".csv", ".tsv", ".txt", ".tab", ".xls", ".xlsx"),
                buttonLabel = "Browse...", placeholder = "No file selected")
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel(
                    "Instructions",
                    h4("Data Upload Instructions"),
                    h5(class="text-info", "Quick Lookup"),
                    p("Enter a UniProt accession number(s) for your protein(s) of interest (e.g. Q01650). Isoform 
                      annotations (e.g. Q01650-1) can be included; however, the specific isoform will not be considered 
                      as SPC scores are indexed by parent protein accession number.  Up to 100 proteins 
                      separated by commas can be searched using this method."),
                    p("If your data are in a form other than UniProt (e.g. ENSEMBL gene, UniGene), a conversion 
                      tool is available", a(href="https://www.uniprot.org/uploadlists/", "here"), "Under 'Select options', 
                      select your ID type in the 'From' field and then 'UniProt KB' in the 'To' field. "),
                    h5(class="text-info", "Bulk Lookup"),
                    p("Upload a csv file containing a single column of UniProt accession numbers, with the header 
                      labeled “Accession”.  Do not include extra characters in the header (e.g. not 'Accession #')."),
                    p("Bulk conversion from a different protein ID type to UniProt is available ",
                      a(href="https://www.uniprot.org/uploadlists/", "here"), 
                      ". Under 'Select options', select your ID type in the 'From' field and then 'UniProt KB'
                      in the 'To' field."),
                    p("With this method, the original upload file will be returned as a downloadable csv file which includes a column containing SPC Scores appended to the original input file.")
                  ),
                  tabPanel(
                    "Quick Lookup",
                    plotOutput("SPC_quick_hist"),
                    tableOutput("SPC_quick_output")
                  ),
                  tabPanel(
                    "Bulk Lookup",
                    plotOutput("SPC_bulk_hist"),
                    tableOutput("SPC_bulk_output"),
                    br(),
                    uiOutput("SPC_csv_dlbutton"),
                    br()
                  )
      )
    )
  ),
  
  ##########    Contact   ##########
  
  tabPanel(
    "Contact",
    div(
      h4(span(class ="text-success", "Contact "), "Us!"),
      p("If you have questions or suggestions for additional features, please contact us by email:"),
      p(class="text-info", style="text-indent:1.5em", "rebekah.gundry at unmc.edu"),
      p("Additional cell surface-related information and tools can be found at our growing website:"),
      p(class="text-info", style="text-indent:1.5em", "www.cellsurfer.net")
    ),
    br()
 ),

  ##########    References   ##########
  
  tabPanel(
    "References",
    div(
      h4("How to reference ", span(class ="text-success", "SurfaceGenie") ),
      p("If you use any of the SurfaceGenie tools in your work, please cite the original manuscript:"),
      p("Waas M, Snarrenberg ST, Littrell J, Jones Lipinski RA, Hansen PA, Corbett JA, Gundry RL, 
        SurfaceGenie: A web-based application for prioritizing cell-type specific marker candidates,", 
        tags$a(href="https://doi.org/10.1101/575969", "https://doi.org/10.1101/575969"))
    ),
    br(),
    div(
      h4("Publications that cite ", span(class ="text-success", "SurfaceGenie") ),
      p("Coming Soon!")
    ),
    br(),
    div(
      h4("Publications that support the ", span(class ="text-success", "SPC Score") ),
      tags$ol(
        tags$li( "Bausch-Fluck D, et al. (2018) The in silico human surfaceome. Proc Natl Acad Sci U S A 115(46):E10988-E10997."),
        tags$li( "da Cunha JP, et al. (2009) Bioinformatics construction of the human cell surfaceome. Proc Natl Acad Sci U S A 106(39):16752-16757"),
        tags$li( "Town J, et al. (2016) Exploring the surfaceome of Ewing sarcoma identifies a new and unique therapeutic target. Proc Natl Acad Sci U S A 113(13):3603-3608" ),
        tags$li( "Diaz-Ramos MC, Engel P, & Bastos R (2011) Towards a comprehensive human cell-surface immunome database. Immunol Lett 134(2):183-187.")
      )
    ),
    br(),
    div(
      h4("Users:"),
      tags$script(type="text/javascript", id="clustrmaps", src="https://cdn.clustrmaps.com/map_v2.js?d=VJztTvZJUQlwpFCwOOYTSK6ktP0YBoNDEMPj1OS_ID0&cl=ffffff&w=a")
      )
    )
  
  ##########    Footer   ##########

#  div(
#    br(), br(),
#    tags$em(p(style="font-size:12px", "Publication Info [Gundry Lab 2018]"))
# )
))
