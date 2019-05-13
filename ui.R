# SurfaceGenie_0.1/ui.R
library(shiny)
library(plotly)

shinyUI(navbarPage("  SurfaceGenie  ", theme = "bootstrap.css",
  
  ##########    Home    ##########
  
  tabPanel(
    "Home",
    div(
     h4("Welcome to ", span(class ="text-success", "SurfaceGenie"),"!"),
     p(tags$i("Integrating predictive and empirical data for rational marker prioritization")),
     p("SurfaceGenie is a tool for analyzing proteomic datasets to identify proteins of 
       interest for immunophenotyping, immunotherapy, drug targeting, and other applications. 
       It works by prioritizing the likelihood that a protein is informative for distinguishing 
       among sample groups (i.e. cell types, experimental conditions). SurfaceGenie generates a 
       score for each protein based on how likely it will be found on the cell surface, the 
       number of samples it is observed in within a comparison set, and the magnitude of the 
       measurement variable (i.e. relative abundance). While a major benefit of SurfaceGenie 
       is the ability to prioritize molecules that are localized to the cell surface, it is 
       also possible to analyze data without this parameter to find proteins of interest that 
       reside in other subcellular localizations. SurfaceGenie works well with approaches that 
       specifically identify cell surface proteins (e.g. Cell Surface Capture) and more generic 
       approaches (e.g. analyses of whole cell lysate). The SurfaceGenie score is context-dependent, 
       meaning that the tool will consider all data within a single dataset input (which may 
       contain multiple experiments and/or cell types). If a user performs a comparison and 
       subsequently determines additional data should be considered, a new file containing 
       all data for the new comparison is required.")
     ),
     p("If you use SurfaceGenie in your research, please cite the article:"),
     p("{ADD Reference and PUBMED LINK HERE}"),
    br(),
    div(
     h4("SurfaceGenie Web Tools"),
     h5(class="text-info", "SurfaceGenie"),
     p(tags$i("Input:  ")),
     p(style="margin-left:1.5em", "SurfaceGenie accepts a .csv file containing a list of proteins 
       (UniProt Accession) and a surrogate value representative of abundance (e.g. number of 
       peptide spectrum matches, peak area) identified within a set of samples. There is no limit 
       to the number of samples that can be analyzed in a single file."),
     p(tags$i("Data Processing:  ")),
     p(style="margin-left:1.5em", "SurfaceGenie calculates the dot product of three independent scores:"),
     tags$ol(
       tags$li(tags$u("Surface Protein Consensus (SPC) score"), br(), "A predictive measure of the likelihood 
               that a particular protein can be present at the cell surface. This value is a sum 
               of the number of predictive datasets for which a protein has been predicted to be 
               localized to the cell surface. Scores range 0-4. For more details on the predictive 
               datasets used, ", a(href="http://google.com", "click here"), "."),
       tags$li(tags$u("Distribution Score"), br(), "A measure of how evenly or unevenly distributed a 
               protein is among multiple samples within a comparison dataset. It is based on the 
               Gini coefficient for calculating statistical dispersion of values. Scores range 0 - 1/(1-N)."),
       tags$li(tags$u("Signal Strength"), br(), "An approximate measure of protein abundance for cell types 
               in which a protein is observed. Proteins at the lower limit of detection are of lower 
               priority than those with more observations, because it is expected that those of higher 
               abundance will practically serve as more accessible markers for downstream technologies. 
               Scores typically range 0 ~ 4 .")
     ),
     h5(class="text-info", "SPC Score Lookup"),
     p("This feature enables users to obtain Surface Protein Consensus (SPC) score for proteins 
       of interest without analyzing data through SurfaceGenie. Users may perform a batch retrieval 
       by uploading a .csv file containing UniProt Accession numbers or may search individual 
       UniProt accession numbers.")
    )
  ),
  
  ##########  SurfaceGenie ##########
  
  tabPanel(
    "SurfaceGenie",
    sidebarPanel(
      h5(class="text-info", "Data Input"),
      fileInput("file1", "Choose CSV File", multiple =FALSE, 
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv"),
                buttonLabel = "Browse...", placeholder = "No file selected"),
      hr(),

      h5(class="text-info", "Scoring Options"),
      #### regular GS
      #### GS, no SPC = iGenie
      #### GS reversed = eineG
      #### GS reversed, no SPC = eineGi
      checkboxGroupInput(
        "scoring_opts", "Select scoring methods:",
        choiceNames = list(
          "SurfaceGenie",
          "eineG",
          "iGenie",
          "eineGi"),
        choiceValues = list(
          "GS", "eineG", "iGenie", "eineGi"),
          selected = list("GS")
      ),
      
      h5(class="text-info", "Processing Options"),
      checkboxGroupInput(
        "processing_opts", "Select processing options:",
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
      hr(),
      h5(class="text-info", "Export Options (for CSV Download Tab)"),
      checkboxGroupInput(
        'export_options', "Choose variables to export:",
        choiceNames = list(
          "SPC score (SPC)",
          "Exclude HLA molecules",
          "CD molecules",
          "Number of CSPA experiments",
          "Gini coefficient (Gini)",
          "Signal strength (SS)",
#          "SurfaceGenie: Genie Score (GS)",
          "UniProt Linkout"),
        choiceValues = list(
#          "SPC", "HLA", "CD", "CSPA #e", "Gini", "SS", "GS", "UniProt Linkout"),
          "SPC", "HLA", "CD", "CSPA #e", "Gini", "SS", "UniProt Linkout"),
        selected = list("GS", "HLA")
      ),
      checkboxInput('bar', 'All/None')
    ),
    
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel(
                    "Instructions",
                    div(
                      h4("Data Upload Instructions"),
                      h5(class="text-info", "Data Format"),
                      p("A .csv file containing a list of proteins (UniProt Accession) and a surrogate value 
                        representative of abundance (e.g. number of peptide spectrum matches, peak area) 
                        identified within a set of samples. "),
                      p("The first column of your data file must be labeled 'Accession' with no extra characters 
                        (e.g. not 'Accession #'). This column should contain the UniProt accession numbers of the 
                        proteins in your samples. You may include isoforms. To convert from a different protein ID 
                        type to UniProt, bulk conversion is available ", a(href="https://www.uniprot.org/uploadlists/", "here"), 
                        ". Under 'Select options', select your ID type in the 'From' field and then 'UniProt KB'
                        in the 'To' field."
                      ),
                      p("Additionally, data files must be in csv format. If you are working in Excel, click 
                        'File --> Save As' and select csv in the drop-down menu to convert from .xlsx to .csv."),
                      h5(class="text-info", "Example Data"),
                      tableOutput("example_data"),
                      br()
                      ),
                    div(
                      h4("Data Proccessing Options"),
                      h5(class="text-info", "Surface Protein Concensus (SPC) Score Consideration"),
                      p("If you are interested in finding cell surface markers, you will want to consider SPC score 
                        when calculating the Genie Score. This is the default setting. If you wish to ignore the SPC 
                        score for your proteins when generating Genie Scores, you may uncheck this option and the SPC 
                        score will be set to 1 for all proteins and will not be weighed into the Genie Score. You can 
                        confirm this in the 'CSV' tab and then uncheck 'SPC' in the export options to remove this from 
                        the download file. This is a feature designed to enable identification of molecules that may 
                        differ among cell types but that may be localized inside the cell."),
                      h5(class="text-info", "HLA Molecule Exclusion"),
                      p("Human leukocyte antigen (HLA) molecules are typically found on the cell surface of most cell 
                        types and due to high sequence similarity among these proteins (e.g. HLA-A3 vs. HLA-A30), it 
                        is often challenging to be certain of the specific gene product based solely on peptide-level 
                        evidence. As a result, it may be useful to exclude these from consideration when attempting to 
                        identify cell surface makers for a specific cell type."),
                      h5(class="text-info", "Find Markers for a Specific Sample"),
                      p("If you are interested in identifying markers that are present in a specific sample (e.g. 
                        positive selection marker for a cell type or experimental condition), SurfaceGenie can exclude 
                        proteins that are not observed in that sample. To do this, select the option “Find markers for 
                        specific sample”. A text box will then appear. In the text box, enter sample name of interest 
                        and make sure it exactly matches what is contained in the file header (i.e. 'd00' for the example 
                        dataset). If you have also selected to have SurfaceGenie group your samples (see 'Sample Grouping 
                        below') then you may also indicate a group (i.e. 'Group 1')."),
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
                      h5(class="text-info", "Plots"),
                      p("Several visualizations are made available by SurfaceGenie:"),
                      tags$ul(
                        tags$li("SurfaceGenie Plot: SurfaceGenie scores plotted in order of priority for all proteins in a dataset."),
                        tags$li("SPC Histogram: Shows the distribution of SPC scores."),
                        tags$li("Clustered Heatmap: Visualize the relationship among samples within a dataset based on the relative 
                                abundance measurement contained in the .csv file."),
                        tags$li("Distribution Score: Shows the distribution of Genie Scores.")
                      ),
                      h5(class="text-info", "CSV File"),
                      p("You may select data to export as columns appended to the right of your original data. 
                        The following variables are available for export:"),
                      tags$ul(
                        tags$li("Surface Protein Concensus score (SPC): A predictive measure of the likelihood that a 
                                particular protein can be present at the cell surface."),
                        tags$li("Distribution Score (Gini): A measure of the distribution of the protein amongst samples. 
                                A higher value corresponds to a more localized distribution.",
                                a(href="https://en.wikipedia.org/wiki/Gini_coefficient", "Wikipedia - Gini coefficient")),
                        tags$li("Signal strength (SS): A weighted value of the maximum value reported among samples for the protein"),
                        tags$li("Genie Score (GS): SurfaceGenie's measure for the value of a protein as a potential 
                                marker of interest."),
                        tags$li("CD molecules (CD): Cluster of differentiation (CD) molecules."),
                        tags$li("Number of CSPA experiments (CSPA-NE): ---"),
                        tags$li("UniProt Linkout: Link to the UniProt for information on the protein.")
                        )
                    )
                  ),
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
                    uiOutput("csv_dlbutton")
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
                    plotlyOutput("SG_dist"),
                    div(class="bnav",
                    uiOutput("SG_dist_PNGdlbutton", class="download_this"),
                    uiOutput("SG_dist_SVGdlbutton", class="download_this")),
                    p(),
                    br(),
                    br(),
                    br(),
                    #eineG
                    plotlyOutput("eineG_dist"),
                    div(class="bnav",
                        uiOutput("eineG_dist_PNGdlbutton", class="download_this"),
                        uiOutput("eineG_dist_SVGdlbutton", class="download_this")),
                    p(),
                    br(),
                    br(),
                    br(),
                    #iGenie
                    plotlyOutput("iGenie_dist"),
                    div(class="bnav",
                        uiOutput("iGenie_dist_PNGdlbutton", class="download_this"),
                        uiOutput("iGenie_dist_SVGdlbutton", class="download_this")),
                    p(),
                    br(),
                    br(),
                    br(),
                    #eineGi
                    plotlyOutput("eineGi_dist"),
                    div(class="bnav",
                        uiOutput("eineGi_dist_PNGdlbutton", class="download_this"),
                        uiOutput("eineGi_dist_SVGdlbutton", class="download_this")),
                    br()
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
      fileInput("file2", "Choose CSV File", multiple =FALSE, 
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv"),
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
                      as SPC scores are indexed by parent protein accession number."),
                    p("If your data are in a form other than UniProt (e.g. ENSEMBL gene, UniGene), a conversion 
                      tool is available", a(href="https://www.uniprot.org/uploadlists/", "here"), "Under 'Select options', 
                      select your ID type in the 'From' field and then 'UniProt KB' in the 'To' field. Up to 100 proteins 
                      separated by commas can be searched using this method."),
                    h5(class="text-info", "Bulk Lookup"),
                    p("Upload a csv file containing a single column of UniProt accession numbers, with the header 
                      labeled “Accession”.  Do not include extra characters in the header (e.g. not 'Accession #')."),
                    p("Bulk conversion from a different protein ID type to UniProt is available ",
                      a(href="https://www.uniprot.org/uploadlists/", "here"), 
                      ". Under 'Select options', select your ID type in the 'From' field and then 'UniProt KB'
                      in the 'To' field."),
                    p("With this method your upload file will be returned as a file available for download which
                      includes a column for SPC scores.")
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
      p(class="text-info", style="text-indent:1.5em", "rgundry at mcw.edu"),
      p("Additional cell surface-related information and tools can be found at our growing website:"),
      p(class="text-info", style="text-indent:1.5em", "www.cellsurfer.net")
    ),
    br()
#    tags$iframe(id = "googleform",
#                src = "https://docs.google.com/forms/d/e/1FAIpQLScRtGpbasA6zokudzm5ujmdatP1bk2AZf_eVIloN8JAftUTVQ/viewform?embedded=true",
#                width = 700,
#                height = 700,
#                frameborder = 0,
#                marginheight = 0)
 ),

  ##########    References   ##########
  
  tabPanel(
    "References",
    div(
      h4("How to reference ", span(class ="text-success", "SurfaceGenie") ),
      p("Latest Version"),
      p(class="text-info", style="text-indent:1.5em", "www.cellsurfer.net"),
      p("Previous Versions"),
      p(class="text-info", style="text-indent:1.5em", "www.cellsurfer.net")
    ),
    br(),
    div(
      h4("Papers referencing ", span(class ="text-success", "SurfaceGenie") ),
      p("Other Cites Go Here"),
      p(class="text-info", style="text-indent:1.5em", "www.cellsurfer.net")
    )
  ),

  ##########    Footer   ##########

  div(
    br(), br(),
    tags$em(p(style="font-size:12px", "Publication Info [Gundry Lab 2018]"))
  )
))
