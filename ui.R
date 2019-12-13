# SurfaceGenie_0.1/ui.R
library(shiny)
library(plotly)

scores <- c("GenieScore", "IsoGenie", "OmniGenie", "IsoOmniGenie")
images <- c("gs.png", "isg.png", "og.png", "iog.png")


shinyUI(navbarPage("", theme = "bootstrap.css",
  
  ##########    Home    ##########
  
  tabPanel(
    "Surface Genie",
    # JavaScript to make links to the "Contact Us" tab
    # Because the anchors are created dynamically, there is no way to know what to use for the 
    # href part of an <a> tag that wants to link to the "Contacts" tab.  Luckily, anchors keep
    # a value called "data-value" that keeps the name of the tab.  So we use JavaScript in the 
    # browser to find all the a tags, look for which one has the data-value set to "contact" 
    # and grab the dynamically generated anchor and click it.
    tags$head(tags$script(HTML('
      var fakeClick = function(tabName) {
        var dropdownList = document.getElementsByTagName("a");
        for (var i = 0; i < dropdownList.length; i++) {
          var link = dropdownList[i];
          if(link.getAttribute("data-value") == tabName) {
            link.click();
          };
        }
      };
      '))),
    
    tags$img(src="website_homepage.png", width="330px", align="right"),
    h4("Welcome to ", span(class ="text-success", "SurfaceGenie"),"!"),
     p(tags$i("Integrating predictive and empirical data for rational marker prioritization")),
     p("SurfaceGenie is a web app for analyzing omic datasets (e.g. proteomic, transcriptomic) to prioritize candidate cell-type specific markers of interest for immunophenotyping, immunotherapy, drug targeting, and other applications. It works by calculating the likelihood a molecule is informative for distinguishing among sample groups (e.g. cell types, experimental conditions)."),
     p("In developing SurfaceGenie we aimed to create an accessible tool for calculation of GenieScore and GenieScore components from input data. Details regarding these calculations can be found in the ", a(href="https://www.biorxiv.org/content/10.1101/575969v2", "corresponding publication"), " or the ", a(href="UserGuide.docx", "User Guide"), ". For all calculations, users are able to export the calculated values and plots generated from analysis of their input data."),
     p("SurfaceGenie was written in R and the web application was developed using the Shiny library. Source code and all reference lookup tables are publicly available ", a(href="https://github.com/GundryLab/SurfaceGenie", "at GitHub"), "."),
     br(), 
#     div(style="width:20%;display:block;margin-left:auto;margin-right:auto",
     div(style="width:25%;display:block;margin-right:auto",
             tags$script(type="text/javascript", id="clustrmaps", src="https://cdn.clustrmaps.com/map_v2.js?d=VJztTvZJUQlwpFCwOOYTSK6ktP0YBoNDEMPj1OS_ID0&cl=ffffff&w=a")
     )
  ),

  
  ##########  Instructions ##########

  tabPanel(
    "Instructions",
    div(
      p(style="font-size: 17px", tags$i("Before you begin:")),
#      p(em("Before you begin:")),
      p("SurfaceGenie contains two separate,  though related,  tools – ", span(class ="text-success", tags$b("GenieScore Calculator")), " and ", span(class ="text-success", tags$b("SPC Score Lookup")), ".  It is strongly recommended that all users read the ", a(href="UserGuide.docx", "User Guide"), " which has step-by-step tutorials for both the GenieScore Calculator and the SPC Score Lookup tools. The User Guide comprehensively defines all of the features available in the SurfaceGenie web application including some background on the theory and calculations."),
      p(style="font-size: 17px", tags$i("Conversion to Uniprot Accession IDs:")),
      p("SurfaceGenie operates with Uniprot Accession IDs only. Bulk conversion of alternate IDs to Uniprot IDs can be performed using the ‘Retrieve/ID mapping tool’ available on the Uniprot website, found here. Note that conversion between IDs is not always one-to-one. Manual curation of the results from the ID mapping is advisable."),
      p(style="font-size: 17px", tags$i("Species availability:")),
      p("Currently, most functions on SurfaceGenie are available only for human, mouse, and rat data. Calculation of some GenieScore permutations do not require Accession numbers, and will work on any type input data (see User Guide for more information). If you have requests for additional species, please ", a(href="#Contact", "contact us", onclick = "fakeClick('Contact')"), "."),
      p(style="font-size: 17px", tags$i("Example files:")),
      p("Examples of files formatted correctly for the GenieScore Calculator and the SPC Score Lookup tools can be downloaded using the links below. For more information, please refer to the ", a(href="UserGuide.docx", "User Guide"), "." ),
      p(a(href="ExampleDataForSurfaceGenie.csv", "GenieScore Calculator example file")),
      p(a(href="ExampleDataForSPCdownload.csv", "SPC Score Lookup example file")),
      br(),
      p(style="font-size: 17px", tags$i("GenieScore Quick Overview:")),
      tags$img(src="GSC_instructions.png", width="800px"), #, align="right"),
      p(style="font-size: 17px", tags$i("SPC Score Lookup Quick Overview:")),
      tags$img(src="SSL_instructions.png", width="800px") #, align="right")
    )
    ),
  
  ##########  SurfaceGenie ##########
  
  tabPanel(
    "GenieScore Calculator",
    sidebarPanel(
      h5(class="text-info", "Data Input"),
      fileInput("file1", "Choose Input File", multiple =FALSE, 
                accept=c(".csv", ".tsv", ".txt", ".tab", ".xls", ".xlsx"), 
                buttonLabel = "Browse...", placeholder = "No file selected"),

      h5(class="text-info", "Species"),
      radioButtons(
        "species", NULL,
        choices = list(
          "Human",
          "Rat",
          "Mouse",
          "Other/Ignore"),
        selected = list("Human")
      ),      
      
      h1(),
      h5(class="text-info", "Scoring Options"),
      
        conditionalPanel(
        condition = "input.species!='Other/Ignore'",
        checkboxGroupInput(
          "scoring_opts", label=NULL,
          choiceNames = mapply(scores, images, FUN=function(score, imgloc) {
            tagList(
              score,
              tags$img(src=imgloc, width=75)
            )
          }, SIMPLIFY = FALSE, USE.NAMES = FALSE),
        choiceValues = list(
          "GS", "IsoGenie", "OmniGenie", "IsoOmniGenie")
        )
      ),
      conditionalPanel(
        condition = "input.species=='Other/Ignore'",
        checkboxGroupInput(
          "scoring_opts_o", label=NULL,
          choiceNames = mapply(scores[c(3,4)], images[c(3,4)], FUN=function(score, imgloc) {
            tagList(
              score,
              tags$img(src=imgloc, width=75)
            )
          }, SIMPLIFY = FALSE, USE.NAMES = FALSE),
          choiceValues = list(
          "OmniGenie", "IsoOmniGenie")
      )
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
      conditionalPanel(
        condition = "input.species!='Other/Ignore'",
        checkboxGroupInput(
        'export_options1', "SurfaceGenie Components:",
        choiceNames = list(
          "SPC score (SPC)",
          "Gini coefficient (Gini)",
          "Signal strength (SS)"
        ),
        choiceValues = list(
          "SPC", "Gini", "SS")
      )
      ),
      conditionalPanel(
        condition = "input.species=='Other/Ignore'",
        checkboxGroupInput(
          'export_options1o', "SurfaceGenie Components:",
          choiceNames = list(
            "Gini coefficient (Gini)",
            "Signal strength (SS)"
          ),
          choiceValues = list(
            "Gini", "SS")
        )
      ),
      conditionalPanel(
        condition = "input.species=='Human'",
        checkboxGroupInput(
        'export_options2h', "Annotations / Link outs:",
        choiceNames = list(
             "HLA molecules",
             "CD molecules",
             "Gene Name",
             "Number of CSPA experiments",
             "Transmembrane",
             "Subcellular Location",
             "UniProt Linkout"),
        choiceValues = list(
          "HLA", "CD", "geneName", "CSPA..e", "Transmembrane", "CC", "UniProt Linkout")
        )
  ),
  conditionalPanel(
    condition = "input.species=='Rat'",
    checkboxGroupInput(
      'export_options2r', "Annotations / Link outs:",
      choiceNames = list(
        "CD molecules",
        "Gene Name",
        "Transmembrane",
        "Subcellular Location",
        "UniProt Linkout"),
      choiceValues = list(
         "CD", "geneName", "Transmembrane", "CC", "UniProt Linkout")
    )
  ),
  conditionalPanel(
    condition = "input.species=='Mouse'",
    checkboxGroupInput(
      'export_options2m', "Annotations / Link outs:",
      choiceNames = list(
        "CD molecules",
        "Gene Name",
        "Number of CSPA experiments",
        "Transmembrane",
        "Subcellular Location",
        "UniProt Linkout"),
      choiceValues = list(
        "CD", "geneName", "CSPA..e", "Transmembrane", "CC", "UniProt Linkout")
    )
  )
  
    ),
    
    mainPanel(
      span(textOutput("txtWarning"), style="color:red"),
      tabsetPanel(type = "tabs",
                  tabPanel(
                    "Data Preview",
                    tableOutput("data_input"),
                    em(textOutput("input_size"))
                  ),
                  tabPanel(
                    "Plots",
                    conditionalPanel(
                      condition = "input.species!='Other/Ignore'",
                      plotOutput("SG_SPC_hist"),
                      div(class="bnav",
                      uiOutput("SG_SPC_hist_PNGdlbutton", class="download_this"),
                      uiOutput("SG_SPC_hist_SVGdlbutton", class="download_this")),
                      p(),
                      br(),
                      br(),
                      br()
                    ),
                    conditionalPanel(
                      condition = "input.species!='Other/Ignore'&input.scoring_opts.indexOf('GS')>-1",
                      plotlyOutput("SG_dist"),
                      div(class="bnav",
                      uiOutput("SG_dist_PNGdlbutton", class="download_this"),
                      uiOutput("SG_dist_SVGdlbutton", class="download_this")),
                      p(),
                      br(),
                      br(),
                      br()
                    ),
                    #IsoGenie
                    conditionalPanel(
                      condition = "input.species!='Other/Ignore'&input.scoring_opts.indexOf('IsoGenie')>-1",
                      plotlyOutput("IsoGenie_dist"),
                      div(class="bnav",
                      uiOutput("IsoGenie_dist_PNGdlbutton", class="download_this"),
                      uiOutput("IsoGenie_dist_SVGdlbutton", class="download_this")),
                      p(),
                      br(),
                      br(),
                      br()
                    ),
                    #OmniGenie
                    conditionalPanel(
                      condition = "(input.species!='Other/Ignore'&input.scoring_opts.indexOf('OmniGenie')>-1)|(input.species=='Other/Ignore'&input.scoring_opts_o.indexOf('OmniGenie')>-1)",
                      plotlyOutput("OmniGenie_dist"),
                      div(class="bnav",
                      uiOutput("OmniGenie_dist_PNGdlbutton", class="download_this"),
                      uiOutput("OmniGenie_dist_SVGdlbutton", class="download_this")),
                      p(),
                      br(),
                      br(),
                      br()
                    ),
                    #IsoOmniGenie
                    conditionalPanel(
                      condition = "(input.species!='Other/Ignore'&input.scoring_opts.indexOf('IsoOmniGenie')>-1)|(input.species=='Other/Ignore'&input.scoring_opts_o.indexOf('IsoOmniGenie')>-1)",
                      plotlyOutput("IsoOmniGenie_dist"),
                      div(class="bnav",
                      uiOutput("IsoOmniGenie_dist_PNGdlbutton", class="download_this"),
                      uiOutput("IsoOmniGenie_dist_SVGdlbutton", class="download_this")),
                      br()
                    )
                  ),
                  tabPanel(
                    "Download Results", 
                    tableOutput("data_output"),
                    conditionalPanel(
                      condition = "output.data_output != undefined",
                      em(textOutput("output_size"))
                    ),
                    br(),
                    uiOutput("csv_dlbutton", class="download_this" ), uiOutput("tsv_dlbutton", class="download_this"), uiOutput("xlsx_dlbutton", class="download_this")
                  )
                  
      )
    )
  ),
  
  ##########  Surface Protein Concensus (SPC) Score  Lookup ##########
  
  tabPanel(
    "SPC Score Lookup",
    sidebarPanel(
      h5(class="text-info", "Input Option 1"),
      textAreaInput("quicklookup", "Uniprot accession number:", 
                placeholder="Enter accession numbers, each on a new line. For example:
                                                              
                                        A0AVT1-1 
                                        A0FGR8-6 
                                        A1L0T0 
                                        A1X283",
                rows=10),
      br(),
      h5(class="text-info", "Input Option 2"),
      fileInput("file2", "Choose Input File", multiple =FALSE, 
                accept=c(".csv", ".tsv", ".txt", ".tab", ".xls", ".xlsx"),
                buttonLabel = "Browse...", placeholder = "No file selected"),
      h1(),
      h5(class="text-info", "Species"),
      radioButtons(
        "species2", NULL,
        choices = list(
          "Human",
          "Rat",
          "Mouse"),
        selected = list("Human")
      )      
      
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
                    "Output Option 1",
                    plotOutput("SPC_quick_hist"),
                    tableOutput("SPC_quick_output")
                  ),
                  tabPanel(
                    "Output Option 2",
                    plotOutput("SPC_bulk_hist"),
                    tableOutput("SPC_bulk_output"),
                    br(),
                    uiOutput("SPC_csv_dlbutton"),
                    br()
                  )
      )
    )
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
    br()
#    div(
#      h4("Users:"),
#      tags$script(type="text/javascript", id="clustrmaps", src="https://cdn.clustrmaps.com/map_v2.js?d=VJztTvZJUQlwpFCwOOYTSK6ktP0YBoNDEMPj1OS_ID0&cl=ffffff&w=a")
#      )
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
)


  
  ##########    Footer   ##########

#  div(
#    br(), br(),
#    tags$em(p(style="font-size:12px", "Publication Info [Gundry Lab 2018]"))
# )
))
