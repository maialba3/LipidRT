# options(shiny.maxRequestSize=500*1024^3) # permitir archivos grandes (hasta 500 GB?) -- Se tiene que poner esta linea antes de lanzar LipidMSapp() si no excede la memoria

options(shiny.maxRequestSize=500*1024^3)
library(shinyjs)

#============================================================================#
# Data
#============================================================================#
##############################################################################
# Reactive variables
rtdataReactive <- reactiveVal(NULL) # rtdata

targetpredictionsReactive <- reactiveVal(NULL)

dbChoicesReactive <- reactiveVal(NULL) # The possible DB values for each lipid class for filtering and prediction.
dbSelectionsReactive <- reactiveVal(NULL) # The selected DB values in the plot for each lipid class for filtering and prediction.

parametersReactive <- reactiveVal(NULL) # Text of parameters displayed on screen in Run.
choicesNamesReactive <- reactiveVal(NULL) # Lipid classes selected for filtering and prediction

# filedbs reactive
filedbspathReactive <- reactiveVal(NULL)
filedbsReactive <- reactiveVal(NULL)
dbsReactive <- reactiveVal(NULL)

adductsTableReactive <- reactiveVal(NULL)
adductsTablepathReactive <- reactiveVal(NULL)


# Guarda las clases seleccionadas
lipid_classes_combined <- reactiveVal(character(0))
selected_lipid_classes_combined <- reactiveVal(character(0))

# The minRT and maxRT of the manual filtering plots, so that the values can be modified
minRT_filter <- reactiveVal(1)
maxRT_filter <- reactiveVal(100)
minC_filter <- reactiveVal(1)
maxC_filter <- reactiveVal(100)

minRT_filter_combined <- reactiveVal(1)
maxRT_filter_combined <- reactiveVal(100)
minC_filter_combined <- reactiveVal(1)
maxC_filter_combined <- reactiveVal(100)

# The minRT and maxRT of the manual prediction plots, so that the values can be modified
minRT_predict <- reactiveVal(1)
maxRT_predict <- reactiveVal(100)
minC_predict <- reactiveVal(1)
maxC_predict <- reactiveVal(100)

selectedclass <- reactiveVal(NULL)

show_filterscore <- reactiveVal(FALSE) # FALSE = all, TRUE = score < rtscorethr

selectedindices_filter <- reactiveVal(NULL)
selectedindices_prediction <- reactiveVal(NULL)
selectedindices_filter_combined <- reactiveVal(NULL)
selectedindices_prediction_combined <- reactiveVal(NULL)

#============================================================================#
# User interface of the Shiny application
#============================================================================#
##############################################################################

ui <- fluidPage(
  useShinyjs(),
  tagList(
    tags$head(
      # Load Font Awesome for icons
      tags$link(rel = "stylesheet", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.5.0/css/all.min.css"),
      tags$style(HTML("

    body {
      padding-top: 10px !important;
    }

/* Change the container panelr */
.shiny-notification-panel {
  right: 20px !important;     /* Keep right-side spacing */
  left: auto !important;      /* Remove any left-side limitation */
  width: auto !important;
  max-width: 90% !important;  /* Allow it to grow up to 90% of the window width */
}

/* General style (in case no type is defined)  */
.shiny-notification {
  background-color: #fc804a;
  border: 1px solid #e74818;
    font-size: 16px !important;
  color: white;
  font-size: 16px;
  padding: 20px;
  max-width: 600px;
  width: 400px;
  height: auto;
  border-radius: 8px;
  margin-right: 20px;
  right: 120px;
  margin-left: auto;
}

/* GREEN (success - message) */
.shiny-notification-message {
  background-color: #669973;
  border: 1px solid #4e7a5d;
}

.shiny-notification-warning {
  background-color: #f5f5f5 !important;
  border: 1px solid #660000 !important;
  color: black !important;
  font-size: 15px !important;
  font-weight: normal !important;
  padding: 20px !important;
  max-width: 1500px !important;       /* Allow it to grow up to 90% of the window width */
  width: 900px !important;
  max-height: 300px !important;
  overflow-y: auto !important;
  border-radius: 8px !important;
  box-shadow: none !important;

    /* Force full opacity */
  opacity: 1 !important;

  /* Important changes to allow it to expand toward the left */
  margin-right: 500px !important;   /* Do not push it to the right */
  right: 620px !important;
}


/* Spinning circle animated loader */
.loader {
  display: inline-block;
  width: 16px;
  height: 16px;
  border: 3px solid white;
  border-top: 3px solid #e74818; /* top border color */
  border-radius: 50%;
  animation: spin 1s linear infinite;
  margin-left: 10px;
  vertical-align: middle;
}

@keyframes spin {
  0% { transform: rotate(0deg); }
  100% { transform: rotate(360deg); }
}


    #console_output {
      overflow-y: auto; /* Enable vertical scrolling if the content overflows. */
      white-space: pre-wrap; /* Adjust the text to preserve the format */
    }

    .dropstart .dropdown-menu {
      position: fixed !important;
      top: 50% !important;
      left: 50% !important;
      transform: translate(-50%, -50%) !important;
      width: 80vw !important; /* % of the total Shiny window it occupies when opened */
      max-width: 1100px;
      height: 95vh;
      overflow: auto;
      z-index: 1050;
      padding: 20px;
      background-color: white;
      box-shadow: 0px 4px 12px rgba(0, 0, 0, 0.15);
      border-radius: 8px;
      justify-content: center;
    }

        /* Checkbox */
        input[type='checkbox'] {
      accent-color: #669999;
        }
        /* Left margin for checkboxes */
        .checkbox {
         margin-left: 15px;
        }






        /* Radio buttons */
    .radio input[type='radio'] {
      accent-color: #669999;
    }
    /* Left margin for Radio buttons */
    .radio {
     margin-left: 15px;
    }


* Visual style of the box */
    .radio label {
      display: inline-block;
      padding: 10px 20px;
      margin-right: 10px;
      border-radius: 6px;
      background-color: transparent;
      border: 1px solid transparent;
      cursor: pointer;
      transition: background-color 0.3s ease;
    }

    /* Selected style */
    .radio label.selected-radio {
      background-color: #f8f8f8;
      border: none;
      box-shadow: none;
      width: 95%;
      font-weight: bold;
    }



    .progress-bar {
      height: 80px;
      font-size: 10px;
      font-weight: bold;
      background-color: #669999;
      color: white;
    }

        /* Slider bar */
    .irs-bar,
    .irs-bar-edge {
    height: 8px !important;
      background: #669999 !important;
      border-color: #669999 !important;
    }



/* Slider handle (redondel) */
.irs-handle {
  width: 12px !important;
  height: 12px !important;
  top: 22px !important;
}

/* Value labels */
.irs-single,
.irs-from,
.irs-to {
  background: #669999 !important;
  border-color: #669999 !important;
  font-size: 10px !important;
  padding: 1px 4px !important;
}

/* Reduce vertical spacing */
.irs {
  margin-top: 0px !important;
  margin-bottom: 0px !important;
  margin-left: 5px !important;
  margin-right: 5px !important;
}




    /* Disable TabPanel */
      .nav-tabs li.disabled a {
        color: grey !important;
        background-color: #f5f5f5 !important;
        cursor: not-allowed !important;
      }


    /* When the checkbox is disabled, the cursor also changes over the text */
    .disabled-checkbox label {
      cursor: not-allowed !important;
      opacity: 0.6;
    }


    table {
      border-collapse: separate !important;
      border-spacing: 50px 0px;
    }

    td {
      max-width: 350px;
      white-space: nowrap;
      overflow: hidden;
      text-overflow: ellipsis;
      padding: 0 15px;
    }
    td:hover {
      overflow: visible;
      white-space: normal;
      background-color: #f9f9f9;
    }

    td.full-message {
      white-space: normal;
      overflow: visible;
      max-width: none;
      background-color: #f9f9f9;
    }


    .no-hover td {
  overflow: hidden !important;
  white-space: normal !important; /* para permitir saltos de línea si quieres */
  background-color: initial !important;
}
.no-hover td:hover {
  overflow: hidden !important;
  white-space: normal !important;
  background-color: initial !important;
}



      .info-box {
        background-color: #fff9e8;
        color: #856404;
        padding: 15px;
        border-radius: 6px;
        display: flex;
        flex-direction: column;
        gap: 10px;
        width: 100%;
        margin-bottom: 30px;
      }
      .info-text {
        font-size: 14px !important;
        line-height: 1.4 !important;
      }

      .toggle-btn {
        width: 100% !important;
        min-height: 35px !important;
        font-size: 11px !important;
        font-weight: bold !important;
        border-color: #856404 !important;
        border-radius: 6px !important;
        cursor: pointer !important;
        margin-bottom: 10px;

        /* Adjustments to allow text wrapping */
        white-space: normal !important; /* allows line breaks */
        word-wrap: break-word !important; /* breaks long words */
        text-align: center; /* optional, centers text horizontally */
        display: flex; /* for vertical centering */
        align-items: center; /* centers content vertically */
        justify-content: center; /* centers content horizontally */
        padding: 0px; /* inner spacing */
      }

      .btn-false {
        background-color: #dedede !important;
        color: black !important;
      }
      .btn-true {
        background-color: #b38b6d !important;
        color: white !important;
      }


  .shinyjs-disabled {
    cursor: not-allowed !important;
    opacity: 0.6;
  }

/* When filedbs_section is disabled */
#filedbs_section.shinyjs-disabled .radio label.selected-radio {
  background-color: transparent !important;
  border: 1px solid transparent !important;
  box-shadow: none !important;
  font-weight: normal !important;
  width: auto !important;
}

/* Remove the visual styling of the radio button */
#filedbs_section.shinyjs-disabled input[type='radio'] {
  accent-color: #cccccc !important;
}

/* Disabled cursor */
#filedbs_section.shinyjs-disabled .radio label {
  cursor: not-allowed !important;
}

#filedbs_section.blocked-section {
  position: relative;
  cursor: not-allowed !important;
}

#filedbs_section.blocked-section::after {
  content: '';
  position: absolute;
  inset: 0;
  background: rgba(255,255,255,0.4);
  cursor: not-allowed;
}



/* === Base style for selectInput === */

/* Background + border */
.form-control {
  background-color: #ffffff;
  border: 1px solid #d0d0d0;
  border-radius: 6px;
  box-shadow: none;
}

/* On focus (when clicked) */
.form-control:focus {
  border-color: #669999;
  box-shadow: 0 0 0 1px rgba(102, 153, 153, 0.3);
}

/* Dropdown arrow color (browser dependent, best effort) */
select.form-control {
  color: #333;
}

/* Reduce vertical height slightly */
select.form-control {
  padding: 4px 8px;
  height: auto;
}


  /* BOTONES tipo toggle */
  .btn-toggle .btn {
    border: 1px solid #669999;
    background: white;
    color: #2c3e50;
    transition: all 0.2s ease;
    border-radius: 6px !important;
  }

  .btn-toggle .btn:hover {
    background: #e6f2f2;
  }

  .btn-toggle .btn.active {
    background: #669999 !important;
    color: white !important;
    box-shadow: none !important;
  }


  /* SWITCH BASE */
  .switch {
    position: relative;
    display: inline-block;
    width: 50px;
    height: 26px;
  }

  .switch input {
    opacity: 0;
    width: 0;
    height: 0;
  }

  /* slider background */
  .slider {
    position: absolute;
    cursor: pointer;
    top: 0; left: 0; right: 0; bottom: 0;
    background-color: #ccc;
    transition: .3s;
    border-radius: 26px;
  }

  /* circle */
  .slider:before {
    position: absolute;
    content: \"\";
    height: 20px;
    width: 20px;
    left: 3px;
    bottom: 3px;
    background-color: white;
    transition: .3s;
    border-radius: 50%;
  }

  /* ON state */
  .switch input:checked + .slider {
    background-color: #669999;
  }

  .switch input:checked + .slider:before {
    transform: translateX(24px);
  }


  /* =========================
     SWITCH CON LABELS
  ========================= */

  .switch-wrapper {
    display: inline-flex;
    align-items: center;
    gap: 6px;
    font-size: 12px;
  }

  .switch-wrapper .form-group {
    margin-bottom: 0 !important;
  }

  .switch-wrapper input[type='checkbox'] {
    width: 52px;
    height: 26px;
    appearance: none;
    background: #ccc;
    border-radius: 26px;
    position: relative;
    cursor: pointer;
    transition: 0.3s;
    outline: none;
  }

  .switch-wrapper input[type='checkbox']:checked {
    background: #669999;
  }

  .switch-wrapper input[type='checkbox']::before {
    content: '';
    position: absolute;
    width: 20px;
    height: 20px;
    top: 3px;
    left: 3px;
    background: white;
    border-radius: 50%;
    transition: 0.3s;
  }

  .switch-wrapper input[type='checkbox']:checked::before {
    transform: translateX(26px);
  }


/* =========================
   SLIDER PERSONALIZADO
========================= */

/* Barra activa */
.irs-bar, .irs-bar-edge {
  background: #669999 !important;
  border-color: #669999 !important;
}

/* Bola (handle) */
.irs-handle {
  border-color: #669999 !important;
  background: #ffffff !important;
}

/* Valor seleccionado */
.irs-single {
  background: #669999 !important;
  color: white !important;
}

/* Tooltip */
.irs-from, .irs-to {
  background: #a3c2c2 !important;
  color: #2c3e50 !important;
}

/* Min / Max */
.irs-min, .irs-max {
  color: #6c757d !important;
}


/* =========================
   RADIO BUTTONS → BOTONES
========================= */

/* Ocultar radio original */
.radio-btn-group input[type='radio'] {
  display: none;
}

/* Espacio entre botones (ajustado) */
.radio-btn-group .radio-inline {
  margin-right: 4px;
}

/* Botón base */
.radio-btn-group label span {
  display: inline-block;
  background-color: #ccc;
  padding: 6px 14px;
  border-radius: 6px;
  cursor: pointer;
  font-weight: 500;
  transition: all 0.15s ease;
}

/* Hover */
.radio-btn-group label span:hover {
  background-color: #b3cccc;
}

/* Activo (efecto pulsado) */
.radio-btn-group input[type='radio']:checked + span {
  background-color: #669999;
  color: white;

  box-shadow:
    inset 0 2px 4px rgba(0,0,0,0.25),
    0 1px 2px rgba(0,0,0,0.2);

  transform: translateY(1px);
}

/* Quitar margen raro */
.radio-btn-group label {
  margin: 0;
}


/* =========================
   CAJAS (reutilizable)
========================= */

.box-green {
  background-color: #e6f2f2;
  border: 1px solid #a3c2c2;
  border-left: 6px solid #669999;
  padding: 10px 12px;
  border-radius: 6px;
  margin-bottom: 10px;
}

.box-grey {
  background-color: #f1f3f5;
  border: 1px solid #ced4da;
  border-left: 6px solid #6c757d;
  padding: 10px 12px;
  border-radius: 6px;
  margin-bottom: 10px;
}



/* =========================
   BOTÓN PRINCIPAL REAL
========================= */
.bootstrap-select > .dropdown-toggle,
.bootstrap-select > .dropdown-toggle.btn-default {
  background-color: #ffffff !important;
  background-image: none !important;
  border: 1px solid #edf5f5 !important;
  color: #2c3e50 !important;
  box-shadow: none !important;
}

/* =========================
   ESTADOS (ESTO ES LO IMPORTANTE)
========================= */
.bootstrap-select.open > .dropdown-toggle,
.bootstrap-select.show > .dropdown-toggle,
.bootstrap-select .dropdown-toggle:focus,
.bootstrap-select .dropdown-toggle:active {
  background-color: #ffffff !important;
  border: 1px solid #edf5f5 !important;
  box-shadow: none !important;
  outline: none !important;
}

/* =========================
   DROPDOWN
========================= */
.dropdown-menu {
  background-color: #ffffff !important;
  border: 1px solid #edf5f5 !important;
}

/* ITEMS */
.dropdown-menu li a {
  color: #2c3e50 !important;
}

/* HOVER */
.dropdown-menu li a:hover {
  background-color: #f7fdfd !important;
  color: #669999 !important;
}

/* SELECTED */
.dropdown-menu .selected a {
  background-color: #f0fbfb !important;
  color: #669999 !important;
}



  ")),
      tags$script(HTML("
      Shiny.addCustomMessageHandler('toggleTabs', function(message) {
        message.ids.forEach(function(id) {
          var tab = $('a[data-value=\"' + id + '\"]');
          if (message.disable) {
            tab.addClass('disabled-tab');
            tab.on('click.disabled', function(e) {
              e.preventDefault();
              e.stopImmediatePropagation();
              return false;
            });
            tab.parent().addClass('disabled');
          } else {
            tab.removeClass('disabled-tab');
            tab.off('click.disabled');
            tab.parent().removeClass('disabled');
          }
        });
      });


      // Apply the class to the label of the selected radio button
      $(document).on('shiny:inputchanged', function(event) {
        if (event.name === 'data_format') {
          // Solo afecta al grupo de data_format
          var $group = $('input[name=\"data_format\"]').closest('.form-group');
          $group.find('label').removeClass('selected-radio');
          $group.find('input[name=\"data_format\"]:checked').closest('label').addClass('selected-radio');
        }

        if (event.name === 'type_chains') {
          // Solo afecta al grupo de type_chains
          var $group = $('input[name=\"type_chains\"]').closest('.form-group');
          $group.find('label').removeClass('selected-radio');
          $group.find('input[name=\"type_chains\"]:checked').closest('label').addClass('selected-radio');
        }

        if (event.name === 'type_lipidclasses') {
          // Solo afecta al grupo de type_lipidclasses
          var $group = $('input[name=\"type_lipidclasses\"]').closest('.form-group');
          $group.find('label').removeClass('selected-radio');
          $group.find('input[name=\"type_lipidclasses\"]:checked').closest('label').addClass('selected-radio');
        }
      });


      // Apply the class when the app loads
      $(document).ready(function() {
        $('input[name=\"data_format\"]:checked').closest('label').addClass('selected-radio');
          $('input[name=\"type_lipidclasses\"]:checked').closest('label').addClass('selected-radio');
          $('input[name=\"type_chains\"]:checked').closest('label').addClass('selected-radio');
      });


      Shiny.addCustomMessageHandler('scrollTop', function(x) {
        var el = document.getElementById(x.id);
        if(el){ el.scrollTop = 0; }
      });
    "))
    ),

    navbarPage(
      position = "static-top",
      theme = shinythemes::shinytheme("flatly"),
      title = " ",
      windowTitle = "LipidRT",
      id = "inTabset",

      # tab1: Import Data
      tabPanel(
        title = "Data import", value = "tab1",
        div(
          style = "display: flex; align-items: center; margin: 0px 0px 20px 20px;",

          # Button image
          tags$img(
              src = "iconohorizontal_lipidrt.png",
              style = "height: auto; width: auto; margin-right: 10px; height: 50px;"
            ),


          tags$h5(
            HTML("<strong style='color:#2C3E50; font-size:18px;'> RT prediction for LC-MS/MS DDA or DIA Data</strong>"),
            style = "margin-right: 30px;"
          ),

        ),

        sidebarPanel(
          width = 4,

          # --- Selector ---
          fluidRow(
            column(12,
                   tags$div(
                     tags$label("Choose your data import method",
                                style = "font-size: 16px; font-weight: bold;"),
                     radioButtons("data_format", label = NULL,
                                  choices = c("Peaklist File" = "data_file",
                                              "RData File" = "data_rdata",
                                              "R environment" = "data_renv"),
                                  selected = "data_file",
                                  width = "100%")
                   ),
                   style = "margin: 0px 0px 0px 0px;"
            )
          ),


          # --- CSV file ---
          conditionalPanel(
            condition = "input.data_format == 'data_file'",

            fluidRow(
              column(12,
                     id = "input_data_file",
                     div(style = "border: 1px solid #ddd; border-radius: 6px;
                         padding: 15px; margin: 10px 20px; background-color: #f9f9f9;",
                         fluidRow(
                           column(6,
                                  tags$label("Column delimiter",
                                             style = "margin-bottom: 0px; font-weight: bold;"),
                                  selectInput("csv_col_sep", NULL,
                                              choices = c(
                                                "Semicolon (;)" = ";",
                                                "Comma (,)" = ",",
                                                "Pipe (|)" = "\\|",
                                                "Tab (\t)" = "\t",
                                                "Space ( )" = " "),
                                              selected = ",")
                           ),
                           column(6,
                                  div(
                                    tags$label("Sheet",
                                               style = "margin-bottom: -2px; font-weight: bold;"),
                                    div(style = "width: 80px;",
                                        tags$style(HTML("
                                        #sheet_file {
                                        height: 35px !important;
                                        padding: 2px 6px !important;
                                        font-size: 14px;
                                        }")),
                                        numericInput("sheet_file", label = NULL,
                                                     value = 1, min = 1, max = 100, step = 1)
                                    )
                                  )
                           )

                         ),

                         fluidRow(
                           column(12,
                                  tags$label("Choose a file to upload",
                                             style = "margin-top: 5px; margin-bottom: 0px; font-weight: bold;"),
                                  uiOutput("file_upload_ui")
                           )
                         )
                     )
              )
            ),
          ),

          # --- RData object ---
          conditionalPanel(
            condition = "input.data_format == 'data_rdata'",

            fluidRow(
              column(12,
                     id = "input_data_rdata",
                     div(style = "border: 1px solid #ddd; border-radius: 6px;
                         padding: 15px; margin: 10px 20px; background-color: #f9f9f9;",
                         fluidRow(
                           column(12,
                                  tags$label("Choose a RData file to upload",
                                             style = "margin-top: 5px; margin-bottom: 0px; font-weight: bold;"),
                                  uiOutput("rdata_upload_ui")
                           )
                         )
                     )
              )
            ),
          ),

          # --- R environment ---
          conditionalPanel(
            condition = "input.data_format == 'data_renv'",
            fluidRow(
              column(12,
                     id = "input_data_renv",
                     div(style = "border: 1px solid #ddd; border-radius: 6px;
                         padding: 15px; margin: 10px 20px; background-color: #f9f9f9;",
                         tags$label("Import from R environment",
                                    style = "font-weight: bold;"),
                         textInput("name_data_renv", NULL, value = "rtdata", width = "100%"),
                         actionButton("import_data_renv", "Import",
                                      style = "color: #fff; background-color: #95a6a6;
                                   border-color: #95a6a6; width: 100%; margin-top: -10px;")
                     )
              )
            )
          ),

          # Button to delete imported data
          fluidRow(
            column(12,
                   div(
                     id = "clear_button_container",
                     style = "text-align: center; margin-top: 10px; display: none;",
                     actionButton("clear_rtdata",
                                  label = tagList(icon("trash-alt"), "Clear imported data"),
                                  style = "background-color: #b3b3b3; color: white; border: none;
                  padding: 8px 18px; font-size: 14px; border-radius: 5px;")
                   )
            )
          ),

          tags$hr(style = "border-color: white"),

          # Allows adding a file containing lipid classes and their adducts for performing annotations and predictions.
          fluidRow(
            column(12,
                   tags$div(
                     tags$label(
                       "Choose your lipid classes source for analysis",
                       style = "font-size: 16px; font-weight: bold;"
                     ),

                     tags$div(
                       "In this section you can modify lipid classes, adducts,
                       and chains that will be used to propose new candidates
                       during RT prediction. To enable these options, you must
                       first import the data; otherwise, this section will
                       remain locked.",
                       style = "color: #7f8c8d; font-size: 13px;
                       margin-top: 2px; margin-bottom: 20px;"
                     ),

                     tags$div(
                       id = "filedbs_section",
                       radioButtons("type_lipidclasses", label = NULL,
                                    choices = list("Default list" = "lipidclasses_default",
                                                   "Custom list"  = "lipidclasses_custom"),
                                    selected = "lipidclasses_default",
                                    width = "100%")
                     )
                   )

            )
          ),

          # If the default list is selected, nothing appears in this tab (the default list is displayed in the annotation and RT prediction tabs, the same in both, synchronized).

          conditionalPanel(
            condition = "input.type_lipidclasses == 'lipidclasses_custom'", # If the customized list is selected…

            fluidRow(
              column(12,
                     div(style = "border: 1px solid #ddd; border-radius: 6px; padding: 15px; margin: 10px 20px; background-color: #f9f9f9;",

                         fluidRow(
                           column(12,
                                  tags$label(
                                    "1) Adducts Table (data frame)",
                                    style = "font-weight: bold; margin-left: 20px;"
                                  ),

                                  div(
                                    style = "border: 1px solid #ddd; border-radius: 6px;
                                    padding: 15px; margin: 0px 0px 35px 0px;
                                    background-color: #FEFEFE;",

                                    # Default vs Custom
                                    radioButtons("type_adductsTable", label = NULL,
                                                 choices = c("Default list" = "adductsTable_default",
                                                             "Custom list" = "adductsTable_custom"),
                                                 selected = "adductsTable_default",
                                                 width = "100%"),


                                    conditionalPanel(
                                      condition = "input.type_adductsTable == 'adductsTable_custom'",

                                      div(
                                        style = "margin: 0px 0px 10px 30px;",
                                        radioButtons("format_adductsTable", label = NULL,
                                                     choices = c("File" = "adductsTable_format_file",
                                                                 "R environment" = "adductsTable_format_renv"),
                                                     selected = "adductsTable_format_file",
                                                     width = "100%"),

                                        conditionalPanel(
                                          condition = "input.format_adductsTable == 'adductsTable_format_file'",

                                          div(
                                            style = "border: 1px solid #ddd; border-radius: 6px;
                                            padding: 15px; margin: 0px 0px 10px -30px;
                                            background-color: #FEFEFE;",
                                            fluidRow(
                                              column(6,
                                                     tags$label("Column delimiter", style = "font-weight: bold;"),
                                                     selectInput("sep_adductsTable", NULL,
                                                                 choices = c("Semicolon (;)" = ";",
                                                                             "Comma (,)" = ",",
                                                                             "Pipe (|)" = "\\|",
                                                                             "Tab (\\t)" = "\t",
                                                                             "Space ( )" = " "),
                                                                 selected = "," )
                                              ),

                                              # Upload
                                              tags$label(
                                                "Choose a adducts table file to upload",
                                                style = "margin-left: 15px; margin-top: 10px; font-weight: bold;"
                                              ),

                                              tags$div(
                                                style = "margin: 0px 15px;",
                                                uiOutput("fileadductsTable_upload_ui")
                                              ),

                                              # download default custom file
                                              tags$div(
                                                style = "text-align: center; margin-top: -15px;",
                                                downloadButton(
                                                  "downloadadductsTable",
                                                  "Default Custom adductsTable.csv",
                                                  class = "btn-primary",
                                                  style = "color: grey; background-color: #ecf0f1;
                                        border-color: #cccccc; width: 60%; padding: 5px 10px;
                                        font-size: 14px; border-width: 1.6px;"
                                                )
                                              )


                                            )
                                          )
                                        ),

                                        conditionalPanel(
                                          condition = "input.format_adductsTable == 'adductsTable_format_renv'",

                                          div(
                                            style = "border: 1px solid #ddd; border-radius: 6px;
                                padding: 15px; margin: 0px 0px 10px -30px;
                                background-color: #FEFEFE;",
                                            textInput("name_adductsTable_renv", NULL,
                                                      value = "adductsTable", width = "100%"),

                                            actionButton("import_adductsTable_renv", "Import adducts table",
                                                         style = "color: #fff; background-color: #95a6a6;
                                border-color: #95a6a6; width: 100%; margin-top: -10px;")
                                          )
                                        )
                                      )
                                    )

                                  )
                           )
                         ),


                         fluidRow(
                           column(12,
                                  tags$label("2) File dbs",
                                             style = "font-weight: bold; margin-left: 20px;"),

                                  div(
                                    style = "border: 1px solid #ddd; border-radius: 6px;
                                    padding: 15px; background-color: #FEFEFE;",
                                    # Default vs Custom
                                    radioButtons("type_filedbs", label = NULL,
                                                 choices = c("Default list" = "filedbs_default",
                                                             "Custom list" = "filedbs_custom"),
                                                 selected = "filedbs_default",
                                                 width = "100%"),

                                    conditionalPanel(
                                      condition = "input.type_filedbs == 'filedbs_custom'",

                                      div(
                                        style = "margin: 0px 0px 10px 30px;",
                                        radioButtons("format_filedbs", label = NULL,
                                                     choices = c("File" = "filedbs_format_file",
                                                                 "R environment" = "filedbs_format_renv"),
                                                     selected = "filedbs_format_file",
                                                     width = "100%"),

                                        conditionalPanel(
                                          condition = "input.format_filedbs == 'filedbs_format_file'",

                                          div(
                                            style = "border: 1px solid #ddd; border-radius: 6px; padding: 15px;
                                            margin: 0px 0px 10px -30px; background-color: #FEFEFE;",

                                            fluidRow(

                                              column(6,
                                                     tags$label("Column delimiter", style = "font-weight: bold;"),
                                                     selectInput("sep_filedbs", NULL,
                                                                 choices = c("Semicolon (;)" = ";",
                                                                             "Comma (,)" = ",",
                                                                             "Pipe (|)" = "\\|",
                                                                             "Tab (\\t)" = "\t",
                                                                             "Space ( )" = " "),
                                                                 selected = ",")
                                              ),

                                              column(6)
                                            ),

                                            tags$label(
                                              "Choose a filedbs file to upload",
                                              style = "margin-left: 0px; margin-top: 10px; font-weight: bold;"
                                            ),

                                            tags$div(
                                              style = "margin: 0px 0px;",
                                              uiOutput("filedbs_upload_ui")
                                            ),

                                            tags$div(
                                              style = "text-align: center; margin-top: -15px;",
                                              downloadButton(
                                                "downloadfiledbs",
                                                "Default Custom filedbs.csv",
                                                class = "btn-primary",
                                                style = "color: grey; background-color: #ecf0f1;
                                                border-color: #cccccc; width: 60%; padding: 5px 10px;
                                                font-size: 14px; border-width: 1.6px;"
                                              )
                                            )

                                          )
                                        ),

                                        conditionalPanel(
                                          condition = "input.format_filedbs == 'filedbs_format_renv'",

                                          div(
                                            style = "border: 1px solid #ddd; border-radius: 6px;
                                            padding: 15px; margin: 0px 0px 10px -30px;
                                            background-color: #FEFEFE;",
                                            textInput("name_filedbs_renv", NULL,
                                                      value = "filedbs", width = "100%"),

                                            actionButton("import_filedbs_renv", "Import filedbs",
                                                         style = "color: #fff; background-color: #95a6a6;
                                              border-color: #95a6a6; width: 100%; margin-top: -10px;")
                                          )
                                        )
                                      )
                                    )

                                  )


                           )
                         )


                     )
              )
            )

          ),

          # Previous and Next buttons
          tags$hr(style = "border-color: white"),
          fluidRow(
            actionButton("JumpTo2", "Next >", width = "100px"),
            style = "margin: 10px 0px 0px 0px"
          )
        ),

        mainPanel(
          conditionalPanel(
            condition = "output.import_summary != null"
          ),
          tableOutput("import_summary") # Here the table with the imported data information will be displayed
        )
      ),

      # tab5: Filter
      tabPanel(title = "Filtering", value = "tab2", # It always appears, regardless of the data type.

               fluidPage(
                 column(6,

                        div(
                          style = "background-color: #fcf9ed; border: 1px solid #ffeeba;
                           border-left: 6px solid #ffdd78;  padding: 15px; border-radius: 6px;
                           margin-bottom: 15px;",

                          div(
                            style = "background-color: #fff3cd; padding: 8px 12px;
                             border-radius: 4px; margin-bottom: 10px; font-weight: bold;
                             color: #7a5d00;",
                            "If your dataset contains a single annotation per feature,
                            set this option to TRUE to allow alternative plasmalogen
                            annotations to be proposed. If multiple annotations
                            per feature are reported, keep it FALSE."
                          ),
                          div(
                            style = "background-color: #ffffff; border: 1px solid #fff3cf;
                           padding: 10px 12px; border-radius: 6px; margin-bottom: 10px;",
                          fluidRow(
                            column(7,
                                   h5(strong("isomeric_plasmalogens")),
                                   h6(
                                     em("If TRUE, alternative plasmalogen annotations
                                    are proposed for features already annotated as
                                    plasmalogens based on m/z matching."),
                                     style = "color:#6c757d;"
                                   )
                            ),

                            column(5,
                                   div(
                                     style = "display: flex; flex-direction: column; align-items: center;",

                                     div(class = "switch-wrapper",
                                         tags$input(
                                           type = "checkbox",
                                           id = "isomeric_plasmalogens"
                                         )
                                     ),
                                     uiOutput("switch_isobaricplasmalogens_status")
                                   )
                            )
                          )
                        )
                        ),

                        tags$hr(),

                        div(
                          style = "background-color: #f8fcfc; border: 1px solid #e2eeee;
                          border-left: 6px solid #ABCAC9; padding: 10px 12px;
                          border-radius: 6px; margin-bottom: 10px;",

                          div(
                            style = "background-color: #E6F2ED; padding: 8px 12px;
                             border-radius: 4px; margin-bottom: 10px; font-weight: bold;
                             color: #3F5F5E;",
                            "Customize your parameters for filtering annotations"
                          ),

                        # Dist C
                        tags$style(HTML("



                         .dist-btn {
                         background-color: #ccc;
                         border: none;
                         padding: 6px 14px;
                         border-radius: 6px;
                         font-weight: 500;
                         transition: 0.2s;
                         }

                         .dist-btn:hover {
                         background-color: #b3cccc;
                         }

                         .dist-btn.active {
                         background-color: #669999 !important;
                         color: white !important;
                         }

                         /* ocultar radio */
                         .radio-btn-group input[type='radio'] {
                         display: none;
                         }

                         /* separación más compacta */
                         .radio-btn-group .radio-inline {
                         margin-right: 0px;
                         }

                         /* botón base */
                         .radio-btn-group label span {
                         display: inline-block;
                         background-color: #ccc;
                         padding: 6px 14px;
                         border-radius: 6px;
                         cursor: pointer;
                         font-weight: 500;
                         transition: all 0.15s ease;
                         }

                         /* hover */
                         .radio-btn-group label span:hover {
                         background-color: #b3cccc;
                         }

                         /* ACTIVO (relieve) 🔥 */
                         .radio-btn-group input[type='radio']:checked + span {
                         background-color: #669999;
                         color: white;
                         /* efecto relieve */
                         box-shadow:
                         inset 0 2px 4px rgba(0,0,0,0.25),
                         0 1px 2px rgba(0,0,0,0.2);
                         transform: translateY(1px); /* efecto pulsado */
                         }

                         /* opcional: quitar margen extra */
                         .radio-btn-group label {
                         margin: 0;
                         }
                                         ")),


                        div(
                          style = "background-color: #ffffff; border: 1px solid #edf5f5;
                           padding: 10px 12px; border-radius: 6px; margin-bottom: 10px;",

                          fluidRow(
                            column(7,
                                   h5(strong("distC_filter")),
                                   h6(
                                     em("Minimum difference in number of carbons from
                                    which lipids are considered distant from the
                                    rest of lipids from the same class. 5 by default."),
                                     style = "color:#2c3e50;"
                                   )
                            ),

                            column(5,
                                   sliderInput("distC_filter", label = NULL,
                                               min = 0, max = 25, value = 5, step = 1)
                            )
                          )
                        ),

                        # minC
                        div(
                          style = "background-color: #ffffff; border: 1px solid #edf5f5;
                           padding: 10px 12px; border-radius: 6px; margin-bottom: 10px;",

                          fluidRow(
                            column(7,
                                   h5(strong("minSpecies_filter")),
                                   h6(
                                     em("Maximum number of lipids allowed in
                                            a distant carbon-number cluster to be
                                            considered for removal. 5 by default"),
                                     style = "color:#2c3e50;"
                                   )
                            ),

                            column(5,
                                   sliderInput("minSpecies_filter", label = NULL,
                                               min = 0, max = 25, value = 5, step = 1)
                            )
                          )
                        ),

                        div(
                          style = "background-color: #ffffff; border: 1px solid #edf5f5;
                           padding: 10px 12px; border-radius: 6px; margin-bottom: 10px;",

                          fluidRow(
                            column(
                              7,
                              h5(strong("distribution")),
                              h6(
                                em("Data distribution for RT modelling. It can be 'gamma' or 'normal'. 'gamma' by default."),
                                style = "color:#2c3e50;"
                              )
                            ),

                            column(
                              5,
                              div(

                                style = "display:flex; justify-content:center; align-items: center; margin-top: 12px;",

                                class = "radio-btn-group",
                                radioButtons(
                                  "distribution_filter",
                                  label = NULL,
                                  choices = c("gamma", "normal"),
                                  selected = "gamma",
                                  inline = TRUE
                                )
                              )
                            )

                          )
                        ),

                        # rtscorethr
                        div(
                          style = "background-color: #ffffff; border: 1px solid #edf5f5;
                           padding: 10px 12px; border-radius: 6px; margin-bottom: 10px;",

                          fluidRow(
                            column(7,
                                   h5(strong("rtscorethr")),
                                   h6(
                                     em("Minimum Score_rtpred value from which previously
                                        filtered lipids will be recovered by the model."),
                                     style = "color:#2c3e50;"
                                   )
                            ),

                            column(5,
                                   sliderInput("score_filter", label = NULL,
                                               min = 0, max = 1, value = 0.5, step = 0.05)
                            )
                          )
                        ))
                        ),
                 column(6,

                        div(
                          style = "background-color: #f8fbfc; border: 1px solid #dde9ee;
                           border-left: 6px solid #9dbbc6;  padding: 15px; border-radius: 6px;
                           margin-bottom: 15px;",

                          div(
                            style = "background-color: #edf5f8; padding: 8px 12px;
                             border-radius: 4px; margin-bottom: 10px; font-weight: bold;
                             color: #47606b;",
                            "Retention-time tolerance for grouping lipid features"
                          ),

                          div(
                            style = "background-color: #ffffff; border: 1px solid #edf5f5;
                           padding: 10px 12px; border-radius: 6px; margin-bottom: 10px;",
                          fluidRow(
                            column(7,
                                   h5(strong("rttol_id")),
                                   h6(
                                     em("Maximum retention-time difference (in seconds)
                                        used to group chemically identical lipid observations.
                                        Observations sharing the same lipid class, subclass
                                        and composition (C and DB) are assigned to the same
                                        GroupID when their retention times fall within ± this
                                        value. For example, a tolerance of 5 seconds allows
                                        grouping observations eluting up to 5 seconds earlier
                                        or 5 seconds later."),
                                     style = "color:#2c3e50;"
                                   )
                            ),

                            column(5,
                                   sliderInput("rttol_id", label = NULL,
                                               min = 0, max = 10, value = 5, step = 0.5)
                            )
                          )
                        )
                        ),

                        div(
                          style = "background-color: #f8fbfc; border: 1px solid #dde9ee;
                           border-left: 6px solid #9dbbc6;  padding: 15px; border-radius: 6px;
                           margin-bottom: 15px;",

                          div(
                            style = "background-color: #edf5f8; padding: 8px 12px;
                             border-radius: 4px; margin-bottom: 10px; font-weight: bold;
                             color: #47606b;",
                            "Retention time tolerance for expanded lipid name assignment"
                          ),

                          div(
                            style = "background-color: #ffffff; border: 1px solid #edf5f5;
                           padding: 10px 12px; border-radius: 6px; margin-bottom: 10px;",
                          fluidRow(
                            column(7,
                                   h5(strong("rttol_name")),
                                   h6(
                                     em("Retention-time tolerance (s) used to match chemically
                                        identical lipids. Annotations with the same lipid class,
                                        subclass and composition within this RT window are considered
                                        equivalent, and the highest-confidence lipid name is assigned."),
                                     style = "color:#2c3e50;"
                                   )
                            ),

                            column(5,
                                   sliderInput("rttol_name", label = NULL,
                                               min = 0, max = 5, value = 1, step = 0.5)
                            )
                          )
                        )
                        ),


                        )

               ),

                         tags$hr(),

                         # Previous and next buttons
                         fluidRow(actionButton("GoBackTo1", "< Previous", width = "100px",
                                               style="margin:10px 0px"),
                                  actionButton("JumpTo3", "Next >", width = "100px",
                                               style="margin:10px 0px"))

      ),

      # tab6: Predict RT
      tabPanel(title = "RT prediction", value = "tab3", # It always appears, regardless of the type of data.

                 fluidPage(
                   column(7,

                          div(
                            style = "background-color: #f8fcfc; border: 1px solid #e2eeee;
                          border-left: 6px solid #ABCAC9; padding: 10px 12px;
                          border-radius: 6px; margin-bottom: 10px;",

                            div(
                              style = "background-color: #E6F2ED; padding: 8px 12px;
                             border-radius: 4px; margin-bottom: 10px; font-weight: bold;
                             color: #27605F;",
                              "Customize your parameters for retention time prediction"
                            ),

                            # Dist C
                            tags$style(HTML("



                         .dist-btn {
                         background-color: #ccc;
                         border: none;
                         padding: 6px 14px;
                         border-radius: 6px;
                         font-weight: 500;
                         transition: 0.2s;
                         }

                         .dist-btn:hover {
                         background-color: #b3cccc;
                         }

                         .dist-btn.active {
                         background-color: #669999 !important;
                         color: white !important;
                         }

                         /* ocultar radio */
                         .radio-btn-group input[type='radio'] {
                         display: none;
                         }

                         /* separación más compacta */
                         .radio-btn-group .radio-inline {
                         margin-right: 0px;
                         }

                         /* botón base */
                         .radio-btn-group label span {
                         display: inline-block;
                         background-color: #ccc;
                         padding: 6px 14px;
                         border-radius: 6px;
                         cursor: pointer;
                         font-weight: 500;
                         transition: all 0.15s ease;
                         }

                         /* hover */
                         .radio-btn-group label span:hover {
                         background-color: #b3cccc;
                         }

                         /* ACTIVO (relieve) 🔥 */
                         .radio-btn-group input[type='radio']:checked + span {
                         background-color: #669999;
                         color: white;
                         /* efecto relieve */
                         box-shadow:
                         inset 0 2px 4px rgba(0,0,0,0.25),
                         0 1px 2px rgba(0,0,0,0.2);
                         transform: translateY(1px); /* efecto pulsado */
                         }

                         /* opcional: quitar margen extra */
                         .radio-btn-group label {
                         margin: 0;
                         }

/* reduce espacio entre opciones */
.checkbox {
  margin-top: 4px;
  margin-bottom: 4px;
}

/* checkbox más limpio */
.checkbox input[type='checkbox'] {
  accent-color: #669999;
}

/* hover suave */
.checkbox label:hover {
  color: #669999;
}
                                         ")),


                            div(
                              style = "background-color: #ffffff; border: 1px solid #edf5f5;
                           padding: 10px 12px; border-radius: 6px; margin-bottom: 10px;",

                              fluidRow(
                                column(7,
                                       h5(strong("dmz")),
                                       h6(
                                         em("Mass tolerance for precursor ions. 5 by default."),
                                         style = "color:#2c3e50;"
                                       )
                                ),

                                # column(5,
                                #        sliderInput("dmz_predict", label = NULL,
                                #                    min = 0, max = 50, value = 5, step = 1)
                                # )
                                column(5,
                                  # uiOutput("dmz_predict_ui")
                                  numericInput(
                                    "dmz_predict",
                                    label = NULL,
                                    value = 5,
                                    min = 0,
                                    step = 1
                                  )
                                )
                              )
                            ),

                            # ppm
                            div(
                              style = "background-color: #ffffff; border: 1px solid #edf5f5;
                           padding: 10px 12px; border-radius: 6px; margin-bottom: 10px;",
                            fluidRow(
                              column(
                                7,
                                h5(strong("ppm")),
                                h6(
                                  em("If TRUE, m/z tolerance in ppm. If FALSE, m/z tolerance in absolute mass error. 'TRUE' by default."),
                                  style = "color:#2c3e50;"
                                )
                              ),

                              column(5,
                                     div(
                                       style = "display: flex; flex-direction: column; align-items: center;",

                                       div(class = "switch-wrapper",
                                           tags$input(
                                             type = "checkbox",
                                             id = "ppm_predict",
                                             checked = "checked"
                                           )
                                       ),
                                       uiOutput("switch_ppm_status")
                                     )
                              )
                            )
                            ),

                            div(
                              style = "background-color: #ffffff; border: 1px solid #edf5f5;
                           padding: 10px 12px; border-radius: 6px; margin-bottom: 10px;",

                              fluidRow(
                                column(7,
                                       h5(strong("distC_predict")),
                                       h6(
                                         em("Minimum difference in number of carbons from
                                    which lipids are considered distant from the
                                    rest of lipids from the same class. 5 by default."),
                                         style = "color:#2c3e50;"
                                       )
                                ),

                                column(5,
                                       sliderInput("distC_predict", label = NULL,
                                                   min = 0, max = 25, value = 5, step = 1)
                                )
                              )
                            ),

                            # minC
                            div(
                              style = "background-color: #ffffff; border: 1px solid #edf5f5;
                           padding: 10px 12px; border-radius: 6px; margin-bottom: 10px;",

                              fluidRow(
                                column(7,
                                       h5(strong("minSpecies_predict")),
                                       h6(
                                         em("Maximum number of lipids allowed in
                                            a distant carbon-number cluster to be
                                            considered for removal. 5 by default"),
                                         style = "color:#2c3e50;"
                                       )
                                ),

                                column(5,
                                       sliderInput("minSpecies_predict", label = NULL,
                                                   min = 0, max = 25, value = 5, step = 1)
                                )
                              )
                            ),

                            # distribution
                            div(
                              style = "background-color: #ffffff; border: 1px solid #edf5f5;
                           padding: 10px 12px; border-radius: 6px; margin-bottom: 10px;",

                              fluidRow(
                                column(
                                  7,
                                  h5(strong("distribution")),
                                  h6(
                                    em("Data distribution for RT modelling. It can be 'gamma' or 'normal'. 'gamma' by default."),
                                    style = "color:#2c3e50;"
                                  )
                                ),

                                column(
                                  5,
                                  div(

                                    style = "display:flex; justify-content:center; align-items: center; margin-top: 12px;",

                                    class = "radio-btn-group",
                                    radioButtons(
                                      "distribution_predict",
                                      label = NULL,
                                      choices = c("gamma", "normal"),
                                      selected = "gamma",
                                      inline = TRUE
                                    )
                                  )
                                )

                              )
                            ),

                            # rtscorethr
                            div(
                              style = "background-color: #ffffff; border: 1px solid #edf5f5;
                           padding: 10px 12px; border-radius: 6px; margin-bottom: 10px;",

                              fluidRow(
                                column(7,
                                       h5(strong("rtscorethr")),
                                       h6(
                                         em("Minimum Score_rtpred value from which a
                                         prediction will be considered valid. 0.5 by default."),
                                         style = "color:#2c3e50;"
                                       )
                                ),

                                column(5,
                                       sliderInput("score_predict", label = NULL,
                                                   min = 0, max = 1, value = 0.5, step = 0.05)
                                )
                              )
                            ),

                            # rt_thr
                            div(
                              style = "background-color: #ffffff; border: 1px solid #edf5f5;
                           padding: 10px 12px; border-radius: 6px; margin-bottom: 10px;",

                              fluidRow(
                                column(7,
                                       h5(strong("rt_thr")),
                                       h6(
                                         em("Threshold for the difference (in seconds) between observed and predicted RT.
                          Applies to new annotations: if the RT score is below the chosen value n rtscorethr and the difference is within this limit,
                          the annotation is marked as discarded but can be manually restored. If it exceeds
                                the limit, it cannot be considered as a new annotation. Default: 30."),
                                         style = "color:#2c3e50;"
                                       )
                                ),

                                column(5,
                                       sliderInput("rt_thr_predict", label = NULL,
                                                   min = 0, max = 60, value = 30, step = 5)
                                )
                              )
                            )

                            )
                   ),

                   # The list of lipid classes and adducts (default or customized)
                   # appears for RT prediction, where you can select which specific
                   # classes you want or not (synchronized with the Annotation list
                   # whenever either of them changes).
                   column(5,
                          div(
                            style = "background-color: #f8fcfc; border: 1px solid #e2eeee;
                          border-left: 6px solid #ABCAC9; padding: 10px 12px;
                          border-radius: 6px; margin-bottom: 10px;",

                            div(
                              style = "background-color: #E6F2ED; padding: 8px 12px;
                             border-radius: 4px; margin-bottom: 10px; font-weight: bold;
                             color: #27605F;",
                              "Selection of lipid subclasses for retention time prediction"
                            ),
                          div(
                            style = "background-color: #ffffff; border: 1px solid #edf5f5;
                           padding: 10px 12px; border-radius: 6px; margin-bottom: 10px;",

                            h5(strong("classes")),

                            h6(
                              em("This section allows selection of lipid subclasses
                                 for retention time (RT) prediction. By default, all
                                 subclasses are selected. Only lipids belonging to
                                 the selected subclasses will be predicted. Use the
                                 dropdown menu to select one or multiple subclasses."),
                              style = "color:#2c3e50;"
                            ),

                            shinyWidgets::pickerInput(
                              "lipidclasses_predict",
                              label = NULL,
                              choices = NULL,
                              multiple = TRUE,
                              options = list(
                                `actions-box` = TRUE,
                                `live-search` = TRUE,
                                `none-selected-text` = "No data selected",
                                `selected-text-format` = "count > 3",
                                `count-selected-text` = "{0} Subclasses selected",
                                `escape` = TRUE
                              )
                            )
                          )
                          )
                   )

                 ),

                 tags$hr(),

                 # Previous and Next buttons
                 fluidRow(
                   actionButton("GoBackTo2", "< Previous", width = "100px", style = "margin:10px 0px"),
                   actionButton("JumpTo4", "Next >", width = "100px", style = "margin:10px 0px")
                 )

      ),

      # tab7: Function Execution
      tabPanel(title = "Run", value = "tab4",

               # LipidMS text and logo at the top
               div(
                 style = "display: flex; align-items: center; margin: 0px 0px 20px 20px;",

                 # Button image
                 tags$img(
                   src = "iconohorizontal_lipidrt.png",
                   style = "height: auto; width: auto; margin-right: 10px; height: 50px;"
                 ),


                 tags$h5(
                   HTML("<strong style='color:#2C3E50; font-size:18px;'> RT prediction for LC-MS/MS DDA or DIA Data</strong>"),
                   style = "margin-right: 30px;"
                 ),

               ),

               # Panel with options to run the analysis
               sidebarPanel(width = 3,
                            fluidRow(
                              # Add the name that will later appear when downloading the files
                              textInput("jobname", "Job Name",
                                        value = paste("Job", as.character(Sys.Date()), sep = "_")),
                              style = "margin: 20px 0px 0px 0px"
                            ),
                            tags$hr(style = "border-color: white"),

                            # # Select the analysis steps to perform (annotation, filtering, and prediction)
                            # fluidRow(h5(HTML("<strong>Select analysis type:</strong>")),
                            #          style = "margin: 20px 0px 0px 0px"),
                            #
                            # fluidRow(
                            #   div(checkboxInput("apply_filter", "Apply Filtering", value = TRUE), style = "margin-left: 20px;"),
                            #   div(checkboxInput("apply_prediction", "Apply RT Prediction", value = TRUE), style = "margin-left: 20px;")
                            # ),
                            #
                            # tags$hr(),

                            # Button to perform the analysis, executes all functions
                            # fluidRow(actionButton("Run_filtering_analysis", "Run Filtering analysis",
                            #                       style="color: #fff; background-color: #669999; border-color: #000; width: 100%;"),
                            #          style = "margin: 20px 0px 0px 0px"),
                            uiOutput("analysis_buttons"),

                            tags$hr(style = "border-color: white"),

                            uiOutput("navButtonsRun")
               ),

               column(9,
                      tags$div(
                        style = "padding: 0px; height: 75vh; overflow-y: auto; position: relative;",

                        # Button to download an .html of the parameters selected so far for the analysis
                        downloadButton("download_parameters_run", "Selected Parameters",
                                       style = "position: sticky; top: 10px; right: 10px; z-index: 900; float: right;"),

                        # Informative text showing the parameters selected so far for the analysis, updating as changes are made
                        uiOutput("parameters_output")
                      )
               )

      ),

      # tab8: Results
      tabPanel(title = "Results", value = "tab5", # If there are no results, this tab is not displayed
               tags$head(
                 tags$style(HTML("
                    .nav-tabs li a {
                      display: block;
                      width: 100%;
                      text-align: center;
                      color: #435f7a;
                      font-weight: bold;
                      margin-bottom: 3px;
                    }

                        .nav-tabs > li.active > a,
                    .nav-tabs > li.active > a:focus,
                    .nav-tabs > li.active > a:hover {
                      color: #2c3e50 !important;
                      background-color: #f9f9f9;
                      margin-bottom: 3px;
                    }

                        .nav-tabs li a {
                      display: block;
                      width: 100%;
                      text-align: center;
                      color: #669999;
                      font-weight: bold;
                    }

                        .nav-tabs > li.active > a,
                    .nav-tabs > li.active > a:focus,
                    .nav-tabs > li.active > a:hover {
                      color: #2c3e50 !important;
                      background-color: #f9f9f9;
                    }


                  #tabs.nav-tabs li a {
                    display: block;
                    width: 100%;
                    text-align: center;
                    color: #435f7a;
                    font-weight: bold;
                    border: 2px solid #f9f9f9 !important;
                    margin-bottom: 3px;
                  }

                  #tabs.nav-tabs > li.active > a,
                  #tabs.nav-tabs > li.active > a:focus,
                  #tabs.nav-tabs > li.active > a:hover {
                    color: #2c3e50 !important;
                    background-color: #f9f9f9;
                    border: 2px solid #2C3E50 !important;
                    margin-bottom: 3px;
                  }


                    .nav-tabs > li > a:hover {
                      color: #2c3e50;
                    }

                    .nav-tabs li {
                      width: 100%;
                    }

                    .no-width-change-tabs .nav-tabs li a {
                      display: block;
                      width: 50%;
                      text-align: center;
                    }
                    .no-width-change-tabs .nav-tabs li {
                      width: 50%;
                    }")),

               ),

               # Sidebar panel with the results tabs
               fluidRow(
                 column(width = 2,

                        # LipidMS text and logo at the top
                        # LipidMS text and logo at the top
                        div(
                          style = "display: flex; align-items: center; margin: 0px 0px 20px 20px;",

                          # tutorial button
                          tags$img(
                            src = "iconohorizontal_lipidrt.png",
                            style = "height: auto; width: auto; margin-right: 10px; height: 50px;"
                          ),

                        ),

                        # Only tabs with data to display will be shown; otherwise, they will not appear
                        sidebarPanel(width = 12,

                                     # SECTION 1: Interactive result exploration
                                     tags$h5(HTML("<strong style='color:#2C3E50; font-size:15px;'> Explore Results</strong>")),
                                     tags$p("Use the tabs above to view peaklists, apply manual filtering, and run predictions.",
                                            style = "color:#7F8C8D; font-size:13px; margin-bottom:10px;"),
                                     # Tabs
                                     tabsetPanel(id = "tabs", type = "tabs",
                                                 tabPanel("Peaklist", value = "peaklist_results"),
                                                 tabPanel("Filter", value = "filter_results"),
                                                 # tabPanel("Filter (Combined View)", value = "filter_combined_results"),
                                                 tabPanel("Prediction", value = "prediction_results"),
                                                 tabPanel("Combined View", value = "filter_combined_results"),
                                     ),

                                     # SECTION 2: Download individual result files
                                     tags$hr(style = "border-color: white"),
                                     tags$h5(HTML("<strong style='color:#2C3E50; font-size:15px;'> Download individual files</strong>")),
                                     fluidRow(
                                       tags$div(uiOutput("download_PeaklistResults_ui"), style = "margin-bottom: 5px;"),
                                       tags$div(uiOutput("download_FilterResults_ui"), style = "margin-bottom: 5px;"),
                                       tags$div(uiOutput("download_PredictionResults_ui"), style = "margin-bottom: 5px;"),
                                       style = "margin: 5px 0px 0px 0px"
                                     ),


                                     # SECTION 3: Download full rtdata object
                                     tags$hr(style = "border-color: white"),
                                     tags$h5(HTML("<strong style='color:#2C3E50; font-size:15px;'> Download full rtdata object</strong>")),
                                     tags$p("Export the entire analysis session, including all results.",
                                            style = "color:#7F8C8D; font-size:13px; margin-bottom:10px;"),
                                     fluidRow(
                                       tags$div(uiOutput("download_rtdata1_ui"), style = "margin-bottom: 5px;"),
                                       tags$div(uiOutput("export_rtdata1_ui")),
                                       style = "margin: 5px 0px 0px 0px"
                                     ),

                                     tags$hr(style = "border-color: white"),

                                     # Navigation buttons
                                     uiOutput("navButtonsResults")
                        )
                 ),

                 # Remaining space for dynamic content that displays the results of the selected tab
                 column(width = 10,
                        conditionalPanel(
                          condition = "input.tabs == 'peaklist_results'", # If 'Batch Results' is selected...
                          uiOutput("features_ui") # displays the feature results
                        ),

                        # Selector LipidClass
                        conditionalPanel(
                          condition = "input.tabs == 'filter_results' || input.tabs == 'prediction_results'",
                          column(
                            width = 6,
                            uiOutput("lipid_class_selector_ui")
                          )
                        ),

                        # Selector LipidClasses (combined view)
                        conditionalPanel(
                          condition = "input.tabs == 'filter_combined_results'",
                          column(
                            width = 6,
                            uiOutput("lipid_classes_combined_selector_ui")
                          )
                        ),

                        conditionalPanel(
                          condition = "input.tabs == 'filter_results'",
                          uiOutput("filter_ui")
                        ),

                        conditionalPanel(
                          condition = "input.tabs == 'filter_combined_results'",
                          uiOutput("filter_combined_ui")
                        ),

                        conditionalPanel(
                          condition = "input.tabs == 'prediction_results'",
                          uiOutput("prediction_ui")
                        )
                 )
               )
      ),

      # tab6: Results
      tabPanel(title = "Target Prediction", value = "tab6", # If there are no results, this tab is not displayed
               tags$head(
                 tags$style(HTML("
                    .nav-tabs li a {
                      display: block;
                      width: 100%;
                      text-align: center;
                      color: #435f7a;
                      font-weight: bold;
                      margin-bottom: 3px;
                    }

                        .nav-tabs > li.active > a,
                    .nav-tabs > li.active > a:focus,
                    .nav-tabs > li.active > a:hover {
                      color: #2c3e50 !important;
                      background-color: #f9f9f9;
                      margin-bottom: 3px;
                    }

                        .nav-tabs li a {
                      display: block;
                      width: 100%;
                      text-align: center;
                      color: #669999;
                      font-weight: bold;
                    }

                        .nav-tabs > li.active > a,
                    .nav-tabs > li.active > a:focus,
                    .nav-tabs > li.active > a:hover {
                      color: #2c3e50 !important;
                      background-color: #f9f9f9;
                    }


                  #tabs.nav-tabs li a {
                    display: block;
                    width: 100%;
                    text-align: center;
                    color: #435f7a;
                    font-weight: bold;
                    border: 2px solid #f9f9f9 !important;
                    margin-bottom: 3px;
                  }

                  #tabs.nav-tabs > li.active > a,
                  #tabs.nav-tabs > li.active > a:focus,
                  #tabs.nav-tabs > li.active > a:hover {
                    color: #2c3e50 !important;
                    background-color: #f9f9f9;
                    border: 2px solid #2C3E50 !important;
                    margin-bottom: 3px;
                  }


                    .nav-tabs > li > a:hover {
                      color: #2c3e50;
                    }

                    .nav-tabs li {
                      width: 100%;
                    }

                    .no-width-change-tabs .nav-tabs li a {
                      display: block;
                      width: 50%;
                      text-align: center;
                    }
                    .no-width-change-tabs .nav-tabs li {
                      width: 50%;
                    }")),

               ),

               # Sidebar panel with the results tabs
               fluidRow(
                 column(width = 2,

                        # LipidMS text and logo at the top
                        div(
                          style = "display: flex; align-items: center; margin: 0px 0px 20px 20px;",

                          # tutorial button
                          tags$img(
                            src = "iconohorizontal_lipidrt.png",
                            style = "height: auto; width: auto; margin-right: 10px; height: 50px;"
                          ),

                        ),

                        # Only tabs with data to display will be shown; otherwise, they will not appear
                        sidebarPanel(width = 12,

                                     # SECTION 1: Interactive result exploration
                                     tags$h5(HTML("<strong style='color:#2C3E50; font-size:15px;'> Explore Results</strong>")),
                                     tags$p("Use the tabs above to view peaklists, apply manual filtering, and run predictions.",
                                            style = "color:#7F8C8D; font-size:13px; margin-bottom:10px;"),
                                     # Tabs
                                     tabsetPanel(id = "tabs1", type = "tabs",
                                                 tabPanel("Targeted Prediction", value = "targeted_results")
                                     ),

                                     # SECTION 2: Download individual result files
                                     tags$hr(style = "border-color: white"),
                                     tags$h5(HTML("<strong style='color:#2C3E50; font-size:15px;'> Download individual files</strong>")),
                                     fluidRow(
                                       tags$div(uiOutput("download_TargetResults_ui"), style = "margin-bottom: 5px;"),
                                       style = "margin: 5px 0px 0px 0px"
                                     ),


                                     # SECTION 3: Download full rtdata object
                                     tags$hr(style = "border-color: white"),
                                     tags$h5(HTML("<strong style='color:#2C3E50; font-size:15px;'> Download full rtdata object</strong>")),
                                     tags$p("Export the entire analysis session, including all results.",
                                            style = "color:#7F8C8D; font-size:13px; margin-bottom:10px;"),
                                     fluidRow(
                                       tags$div(uiOutput("download_rtdata2_ui"), style = "margin-bottom: 5px;"),
                                       tags$div(uiOutput("export_rtdata2_ui")),
                                       style = "margin: 5px 0px 0px 0px"
                                     ),

                                     tags$hr(style = "border-color: white"),

                                     # Navigation buttons
                                     fluidRow(
                                       actionButton("GoBackTo5", "< Previous", width = "100px", style = "margin:10px 0px"),
                                       style = "margin: 0px 0px 0px 0px"
                                     )
                        )
                 ),

                 # Remaining space for dynamic content that displays the results of the selected tab
                 column(width = 10,
                        conditionalPanel(
                          condition = "input.tabs1 == 'targeted_results'", # If 'Batch Results' is selected...
                          uiOutput("targeted_ui") # displays the feature results
                        )
                 )
               )
      )
    )
  )
)



#============================================================================#
# Shiny application server
#============================================================================#
##############################################################################
server <- function(input, output, session) {

  #============================================================================#
  # Button functionality: Next, Previous
  #============================================================================#
  observeEvent(input$JumpTo2, {updateTabsetPanel(session, "inTabset", selected = "tab2")}) # Data Import -> Filtering
  observeEvent(input$GoBackTo1, {updateTabsetPanel(session, "inTabset", selected = "tab1")}) # Filtering -> Data Import
  observeEvent(input$JumpTo3, {updateTabsetPanel(session, "inTabset", selected = "tab3")}) # Filtering -> Prediction
  observeEvent(input$GoBackTo2, {updateTabsetPanel(session, "inTabset", selected = "tab2")}) # Prediction -> Filtering
  observeEvent(input$JumpTo4, {updateTabsetPanel(session, "inTabset", selected = "tab4")}) # Prediction -> Run
  observeEvent(input$GoBackTo3, {updateTabsetPanel(session, "inTabset", selected = "tab3")}) # Run -> Prediction
  observeEvent(input$JumpTo5, {updateTabsetPanel(session, "inTabset", selected = "tab5")}) # Run -> Results
  observeEvent(input$GoBackTo4, {updateTabsetPanel(session, "inTabset", selected = "tab4")}) # Results -> Run
  observeEvent(input$JumpTo6, {updateTabsetPanel(session, "inTabset", selected = "tab6")}) # Results -> Targeted
  observeEvent(input$GoBackTo5, {updateTabsetPanel(session, "inTabset", selected = "tab5")}) # Targeted -> Results


  #============================================================================#
  # Run Tab
  #============================================================================#
  output$analysis_buttons <- renderUI({

    filtering_done <- !is.null(rtdataReactive()$filtering$data)

    tagList(

      fluidRow(
        # actionButton(
        #   "Run_filtering_analysis",
        #   "Run Filtering Analysis",
        #   style = "font-weight: bold; font-size:15px; color: #fff; background-color: #669999; border-color: #000; width: 100%;"
        # ),
        actionButton(
          "Run_filtering_analysis",
          HTML('<span style="
        display:inline-block;
        width:20px;
        height:20px;
        line-height:20px;
        border-radius:50%;
        background-color:white;
        color:#669999;
        text-align:center;
        font-weight:bold;
        margin-right:8px;
      ">1</span> Run Filtering Analysis'),
          style = "font-weight: bold; font-size:15px; color: #fff; background-color: #669999; border-color: #000; width: 100%;"
        ),
        style = "margin: 20px 0px 0px 0px"
      ),

      div(
        style = "
        margin-top:10px;
        padding:10px;
        background-color:#f8f9fa;
        border-left:4px solid #669999;
        font-size:13px;
      ",
        HTML(
          "<b>Filtering step:</b> Clicking <i>Run Filtering Analysis</i> will perform the filtering process.
        Once the analysis is finished, go to <b>Results → Filter</b> and carefully review the filtering results.
        If everything looks correct, return to the <b>Run</b> tab and proceed with the prediction analysis."
        )
      ),

      if (filtering_done) {

        tagList(

          br(),

          fluidRow(
            # actionButton(
            #   "Run_prediction_analysis",
            #   "Run Prediction Analysis",
            #   # style = "color: #fff; background-color: #336699; border-color: #000; width: 100%;"
            #   style = "font-weight: bold; font-size:15px; color: #fff; background-color: #669999; border-color: #000; width: 100%;"
            # ),
            actionButton(
              "Run_prediction_analysis",
              HTML('<span style="
        display:inline-block;
        width:20px;
        height:20px;
        line-height:20px;
        border-radius:50%;
        background-color:white;
        color:#669999;
        text-align:center;
        font-weight:bold;
        margin-right:8px;
      ">2</span> Run Prediction Analysis'),
              style = "font-weight: bold; font-size:15px; color: #fff; background-color: #669999; border-color: #000; width: 100%;"
            ),
            style = "margin: 20px 0px 0px 0px"
          ),

          div(
            style = "
        margin-top:10px;
        padding:10px;
        background-color:#f8f9fa;
        border-left:4px solid #669999;
        font-size:13px;
      ",
            HTML(
              "<b>Prediction step:</b> Filtering results have been detected.
            After confirming that the filtering output is correct, click
            <i>Run Prediction Analysis</i> to generate the RT predictions."
            )
          )
        )
      }
    )
  })


  #============================================================================#
  # Hide Tabs
  #============================================================================#
  observe({
    rtdata <- rtdataReactive()

    if (is.null(rtdata)) {
      hideTab("inTabset", "tab5")
      hideTab("inTabset", "tab6")
      output$navButtonsRun <- renderUI({
        fluidRow(
          actionButton("GoBackTo3", "< Previous", style = "margin:10px 10px"),
          style = "margin: 0px 0px 0px 0px"
        )
      })

    } else {
      showTab("inTabset", "tab5")
      output$navButtonsResults <- renderUI({
        fluidRow(
          actionButton("GoBackTo4", "< Previous", style = "margin:10px 10px"),
          style = "margin: 0px 0px 0px 0px"
        )
      })

      if (is.null(rtdata$filtering$data)) { # if there is no filtering...
        hideTab("tabs", "filter_results")
        hideTab("tabs", "filter_combined_results")

      } else {
        showTab("tabs", "filter_results")
        showTab("tabs", "filter_combined_results")

        output$navButtonsResults <- renderUI({
          fluidRow(
            actionButton("GoBackTo4", "< Previous", style = "margin:10px 10px"),
            actionButton("JumpTo6", "Next >", style = "margin:10px 10px"),
            style = "margin: 0px 0px 0px 0px"
          )
        })
      }
      if (is.null(rtdata$prediction$candidates)) { # if there are no predictions...
        hideTab("tabs", "prediction_results")
        hideTab("inTabset", "tab6") # target
      } else {
        showTab("tabs", "prediction_results")
        showTab("inTabset", "tab6") # target
      }

      output$navButtonsRun <- renderUI({
        fluidRow(
          actionButton("GoBackTo3", "< Previous", style = "margin:10px 10px"),
          actionButton("JumpTo5", "Next >", style = "margin:10px 10px"),
          style = "margin: 0px 0px 0px 0px"
        )
      })
    }
  })



  ##############################################################################
  # Filtering parameters
  ##############################################################################
  output$switch_isobaricplasmalogens_status <- renderUI({
    if (isTRUE(input$isomeric_plasmalogens)) {
      tags$span(
        style = "font-weight: bold; color: #669999;",
        "TRUE"
      )
    } else {
      tags$span(
        style = "font-weight: bold; color: black;",
        "FALSE"
      )
    }
  })

  ##############################################################################
  # Filtering parameters
  ##############################################################################
  output$switch_ppm_status <- renderUI({
    if (isTRUE(input$ppm_predict)) {
      tags$span(
        style = "font-weight: bold; color: #669999;",
        "TRUE"
      )
    } else {
      tags$span(
        style = "font-weight: bold; color: black;",
        "FALSE"
      )
    }
  })



  #============================================================================#
  # Upload file Outputs "Data Import"
  #============================================================================#
  output$file_upload_ui <- renderUI({
    fileInput("path_data_file", NULL,
              multiple = FALSE,
              accept = c("text/csv", "text/plain", ".csv", ".txt", ".xlsx"),
              width = "100%",
              placeholder = "No file selected")
  })

  output$rdata_upload_ui <- renderUI({
    fileInput("path_data_rdata", NULL,
              multiple = FALSE,
              accept = c("text/plain", ".RData", ".gz"),
              width = "100%",
              placeholder = "No file selected")
  })

  output$fileadductsTable_upload_ui <- renderUI({
    fileInput("fileadductsTable", NULL,
              multiple = FALSE,
              accept = c("text/plain", ".csv"),
              width = "100%",
              placeholder = "No file selected")
  })

  output$filedbs_upload_ui <- renderUI({
    fileInput("filedbs", NULL,
              multiple = FALSE,
              accept = c("text/plain", ".csv"),
              width = "100%",
              placeholder = "No file selected")
  })





  observe({
    if (is.null(rtdataReactive())) {
      shinyjs::disable("filedbs_section")
      shinyjs::addClass("filedbs_section", "blocked-section")
    } else {
      shinyjs::enable("filedbs_section")
      shinyjs::removeClass("filedbs_section", "blocked-section")
    }
  })


    # -----------------------------
    # Block / unblock inputs
    # -----------------------------
    observe({
      if (!is.null(rtdataReactive())) {
        shinyjs::show("clear_button_container")

        # Block everything except the delete button
        shinyjs::disable("data_format")
        shinyjs::disable("input_data_file")
        shinyjs::disable("input_data_rdata")
        shinyjs::disable("input_data_renv")
      } else {
        shinyjs::hide("clear_button_container")

        # Unblock everything
        shinyjs::enable("data_format")
        shinyjs::enable("input_data_file")
        shinyjs::enable("input_data_rdata")
        shinyjs::enable("input_data_renv")
      }
    })

    # ------------------------------------
    # Button to clear imported data
    # ------------------------------------
    observeEvent(input$clear_rtdata, {
      rtdataReactive(NULL)
      showNotification("Imported data has been cleared.", type = "message")

      if (input$data_format == "data_file") {
        if (!is.null(input$path_data_file)) {
          output$file_upload_ui <- renderUI({
            fileInput("path_data_file", NULL,
                      multiple = FALSE,
                      accept = c("text/csv", "text/plain", ".csv", ".txt", ".xlsx"),
                      width = "100%",
                      placeholder = "No file selected")
          })
        }

      } else if (input$data_format == "data_rdata") {
        if (!is.null(input$path_data_rdata)) {
          output$rdata_upload_ui <- renderUI({
            fileInput("path_data_rdata", NULL,
                      multiple = FALSE,
                      accept = c("text/plain", ".RData", ".gz"),
                      width = "100%",
                      placeholder = "No file selected")
          })
        }

      } else if (input$data_format == "data_renv") {
        if (input$name_data_renv != "rtdata") {
          updateTextInput(session, "name_data_renv", value = "rtdata")
        }
      }
    })



  #============================================================================#
  # Show lipids score < threshold
  #============================================================================#
  output$toggle_button_ui <- renderUI({
    if (show_filterscore()) {
      actionButton("toggle_btn", "Show all lipids", class = "toggle-btn btn-true")
    } else {
      actionButton("toggle_btn", "Show lipids with score < threshold", class = "toggle-btn btn-false")
    }
  })

  observeEvent(input$toggle_btn, {
    show_filterscore(!show_filterscore())
  })


  #============================================================================#
  # Dynamic text: Selected parameters from Run
  #============================================================================#
  output$parameters_output <- renderUI({

    # tittle: Selected parameters
    text_output <- paste0(
      "<h4 style='color:#669999; font-weight:bold; margin-bottom:15px'>",
      "<i class='fa-solid fa-circle-info' style='color:#669999; margin-right:8px;'></i>",
      "Selected Parameters</h4>",

      "<p style='color:#666666; font-style:italic; margin-top:25px; margin-bottom:25px;'>",
      "If you want to download the selected parameters as an <b>.html</b> file, make sure to do so <b>before clicking the 'Run'</b> button.<br>",
      "Otherwise, the displayed values may change dynamically during the analysis.",
      "</p>"
    )


    if (!is.null(rtdataReactive())) {
      # Data source:
      source <- rtdataReactive()$metadata$source
      polarity <- rtdataReactive()$metadata$polarity
    }


    # Choices: Polarity
    choices_polarity <- c(
      "positive" = "positive",
      "negative" = "negative"
    )

text_output <- paste0(
  text_output,
  "<table style='border-collapse: collapse; margin-left: 10px; margin-bottom: 10px;'>",

  # Job name
  "<tr>",
  "<td style='font-weight:bold; color:#2c3e50; padding-right:30px;'>●&nbsp;&nbsp; Job name</td>",
  "<td style='color:#2c3e50;'>", input$jobname, "</td>",
  "</tr>",

  # If no data has been imported
  if (is.null(rtdataReactive())) paste0(
    "<tr>",
    "<td colspan='2' class='full-message' style='color:#B03A2E; font-style:italic; padding-top:10px;'>⚠️ No data has been uploaded or imported yet. <br>Please ensure you have imported a valid file or environment object.</td>",
    "</tr>"
  ) else paste0(
    # File/Object name depending on data_format
    if (input$data_format == "data_file") paste0(
      "<tr>",
      "<td style='font-weight:bold; color:#2c3e50; padding-right:30px;'>●&nbsp;&nbsp; File name</td>",
      "<td style='color:#2c3e50;'>", input$path_data_file$name, "</td>",
      "</tr>",
      "<tr>",
      "<td style='font-weight:bold; color:#2c3e50; padding-right:30px;'>●&nbsp;&nbsp; File sheet</td>",
      "<td style='color:#2c3e50;'>", input$sheet_file, "</td>",
      "</tr>"
    ) else if (input$data_format == "data_renv") paste0(
      "<tr>",
      "<td style='font-weight:bold; color:#2c3e50; padding-right:30px;'>●&nbsp;&nbsp; Object name</td>",
      "<td style='color:#2c3e50;'>", input$name_data_renv, "</td>",
      "</tr>"
    ),

    # Data source
    "<tr>",
    "<td style='font-weight:bold; color:#2c3e50; padding-right:30px;'>●&nbsp;&nbsp; Data source</td>",
    "<td style='color:#2c3e50;'>", source, "</td>",
    "</tr>",

    # Polarity
    "<tr>",
    "<td style='font-weight:bold; color:#2c3e50; padding-right:30px;'>●&nbsp;&nbsp; Polarity</td>",
    "<td style='color:#2c3e50;'>", names(choices_polarity)[choices_polarity == polarity], "</td>",
    "</tr>"
  ),

  "</table>",
  "<hr style='border:none; height:1px; background-color:#ffffff;'>"
)

    # Choices: type list lipid classes (default or custom)
    if (input$type_lipidclasses == "lipidclasses_default") {
      text_type_lipidclasses<- "&nbsp;&nbsp; (default list)"
    } else if (input$type_lipidclasses == "lipidclasses_custom") {
      text_type_lipidclasses <- "&nbsp;&nbsp; (custom list)"
    }


    # --- Filtering ---
    text_output <- paste0(text_output,
                          "<p><hr style='border:none; height:1px; background-color:#ffffff;'></p>")

    # Filtering (TRUE)
    # if (input$apply_filter) {
    if (TRUE) {
      text_output <- paste0(text_output,
                            "<p><i class='fa-solid fa-square-check' style='color:#95a5a6;'></i> <b> Filtering</b></p>",

                            "<table style='border-collapse: collapse; margin-left: 50px;'>",

                            "<tr>",
                            "<th style='padding-right: 50px;'></th>",
                            "<th'></th>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; distC_filter</td>",
                            "<td>", input$distC_filter, " C</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; minSpecies_filter</td>",
                            "<td>", input$minSpecies_filter, " C</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; distribution</td>",
                            "<td>", input$distribution_filter, "</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; score</td>",
                            "<td>", input$score_filter, "</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; isomeric_plasmalogens</td>",
                            "<td>", input$isomeric_plasmalogens, "</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; rttol_id</td>",
                            "<td> ± ", input$rttol_id, " s</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; rttol_name</td>",
                            "<td> ± ", input$rttol_name, " s</td>",
                            "</tr>",

                            "</table>"

      )
    } else { # Filtering (FALSE)
      text_output <- paste0(text_output,
                     "<p><i class='fa-regular fa-square' style='color:gray;'></i> <b> Filtering</b></p>",
                     "<p style='margin-left:20px; color:gray;'>Filtering will not be applied.</p>"
      )
    }


    # --- Prediction ---
    text_output <- paste0(text_output,
                          "<p><hr style='border:none; height:1px; background-color:#ffffff;'></p>")
    # Prediction (TRUE)
    # if (input$apply_prediction) {
    if (TRUE) {
      text_output <- paste0(text_output,
                            "<p><i class='fa-solid fa-square-check' style='color:#95a5a6;'></i> <b> RT Prediction</b></p>",

                            "<table style='border-collapse: collapse; margin-left: 50px;'>",

                            "<tr>",
                            "<th style='padding-right: 50px;'></th>",
                            "<th'></th>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; dmz</td>",
                            "<td>", input$dmz_predict, "</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; ppm</td>",
                            "<td>", input$ppm_predict, "</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; distC_predict</td>",
                            "<td>", input$distC_predict, " C</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; minSpecies_predict</td>",
                            "<td>", input$minSpecies_predict, " C</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; distribution</td>",
                            "<td>", input$distribution_predict, "</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; score</td>",
                            "<td>", input$score_predict, "</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; rt_thr</td>",
                            "<td>", input$rt_thr_predict, " s</td>",
                            "</tr>",

                            "</table>"

      )

      # dropdown row with the lipid classes to predict selected
      if (!is.null(filedbsReactive())) {
        filedbs <- filedbsReactive()


        if (!is.null(rtdataReactive())) {

          if (rtdataReactive()$metadata$polarity == "positive") {
            show_pos <- TRUE
            show_neg <- FALSE
          } else if (rtdataReactive()$metadata$polarity == "negative") {
            show_pos <- FALSE
            show_neg <- TRUE
          } else {
            show_pos <- TRUE
            show_neg <- TRUE
          }

        } else {
          show_pos <- TRUE
          show_neg <- TRUE
        }




        lipidclasses <- unique(filedbs$classes)
        # selected lipid classes by the user
        selected_lipids <- lipidclasses[lipidclasses %in% input$lipidclasses_predict]

        unique_id <- paste0("lipid_details_", as.integer(runif(1, 1000, 9999)))  # unique ID

        if (length(selected_lipids) == 0) {
          lipid_text <- "<div style='margin-left: 30px; color:gray; font-style: italic;'>No lipid classes selected for RT prediction.</div>"
        } else {
          lipid_text <- "<table style='border-collapse: collapse; margin-left: 30px; margin-top:10px;'>"


          lipid_text <- paste0(
            lipid_text,
            "<table style='border-collapse: collapse; border-spacing: 10px 5px; font-weight:bold; color:#2c3e50;'>",
            "<thead><tr>",
            "<th style='padding: 1px 15px 1px 0;'>classes</th>",
            if (show_pos) "<th style='padding: 1px 0;'>adducts_pos</th>" else "",
            if (show_neg) "<th style='padding: 1px 0;'>adducts_neg</th>" else "",
            "<th style='padding: 1px 0;'>fas</th>",
            "<th style='padding: 1px 0;'>sph</th>",
            "<th style='padding: 1px 0;'>alkyl</th>",
            "<th style='padding: 1px 0;'>vinyl</th>",
            "</tr></thead><tbody>"
          )


          for (lipid in selected_lipids) {

            if (show_pos) {
              adductspos_text <- filedbs$adducts_pos[filedbs$classes == lipid]
              adductspos_text <- gsub("([;,|])\\s*", "\\1 ", adductspos_text)
            }

            if (show_neg) {
              adductsneg_text <- filedbs$adducts_neg[filedbs$classes == lipid]
              adductsneg_text <- gsub("([;,|])\\s*", "\\1 ", adductsneg_text)
            }


            fas_text <- filedbs$fas[filedbs$classes == lipid]
            fas_text <- gsub("([;,|])\\s*", "\\1 ", fas_text)

            sph_text <- filedbs$sph[filedbs$classes == lipid]
            sph_text <- gsub("([;,|])\\s*", "\\1 ", sph_text)

            alkyl_text <- filedbs$alkyl[filedbs$classes == lipid]
            alkyl_text <- gsub("([;,|])\\s*", "\\1 ", alkyl_text)

            vinyl_text <- filedbs$vinyl[filedbs$classes == lipid]
            vinyl_text <- gsub("([;,|])\\s*", "\\1 ", vinyl_text)

            # Compact row with lipid in bold
            lipid_text <- paste0(
              lipid_text,
              "<tr>",
              "<td style='padding: 1px 15px 1px 0; font-size: 85%;'><b>", lipid, "</b></td>",
              if (show_pos) paste0("<td style='padding: 1px 0; font-size: 85%;'>", adductspos_text, "</td>") else "",
              if (show_neg) paste0("<td style='padding: 1px 0; font-size: 85%;'>", adductsneg_text, "</td>") else "",
              "<td style='padding: 1px 0; font-size: 85%;'>", fas_text, "</td>",
              "<td style='padding: 1px 0; font-size: 85%;'>", sph_text, "</td>",
              "<td style='padding: 1px 0; font-size: 85%;'>", alkyl_text, "</td>",
              "<td style='padding: 1px 0; font-size: 85%;'>", vinyl_text, "</td>",
              "</tr>"
            )

          }

          lipid_text <- paste0(lipid_text, "</table>")
        }

        # Expandable
        expandable_row <- paste0(
          "<div style='margin-left: 50px; margin-top: 20px;'>",
          "<div style='color:#2c3e50; font-weight:bold; cursor:pointer;' onclick='toggleDetails(\"", unique_id, "\")'>",
          "<span id='", unique_id, "_toggle'>△</span> Show lipid classes to predict",
          "<span style='font-weight:normal;'>", text_type_lipidclasses, "</span></div>",
          "<div id='", unique_id, "' style='display:block; margin-top: 10px; color:#444;'>",
          lipid_text,
          "</div>",
          "</div>"
        )

        # Output final
        text_output <- paste0(
          text_output,
          expandable_row,
          "<script>
      function toggleDetails(id) {
        var content = document.getElementById(id);
        var toggleIcon = document.getElementById(id + '_toggle');
        if (content.style.display === 'none') {
          content.style.display = 'block';
          toggleIcon.innerHTML = '△';
        } else {
          content.style.display = 'none';
          toggleIcon.innerHTML = '▽';
        }
      }
    </script>"
        )
      }

    } else { # Prediction (FALSE)
      text_output <- paste0(text_output,
                            "<p><i class='fa-regular fa-square' style='color:gray;'></i> <b> Prediction</b></p>",
                            "<p style='margin-left:20px; color:gray;'>RT Prediction will not be applied.</p>"
      )
    }

    parametersReactive(text_output)
    HTML(text_output)
  })


  #============================================================================#
  # Download .html with the selected parameters
  #============================================================================#
  output$download_parameters_run <- downloadHandler(
    filename = function(){paste(input$jobname, "parameters.html", sep="_")},
    content = function(file) {
      raw_html <- parametersReactive()

      full_html <- paste0(
        "<!DOCTYPE html>\n<html>\n<head>\n",
        "<meta charset='UTF-8'>\n",
        "<link rel='stylesheet' href='https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css'>\n",
        "<link rel='stylesheet' href='https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.5.0/css/all.min.css'>\n",
        "<style>\n",
        "body { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; font-size: 14px; color: #2c3e50; padding: 30px; margin-left: 50px; }\n",

        "table {border-collapse: separate !important; /* importante */
        border-spacing: 50px 0px;  /* 20px horizontal, 5px vertical */
        }\n",
        "td {max-width: 350px; white-space: nowrap; overflow: hidden;
          text-overflow: ellipsis; adding: 0 15px; /* 0 arriba/abajo, 15px izquierda/derecha */
        }\n",
        "td:hover {overflow: visible; white-space: normal; background-color: #f9f9f9;
        }\n",

        "</style>\n",
        "</head>\n<body>\n",
        raw_html,
        "\n</body>\n</html>"
      )

      writeLines(full_html, file)
    }
  )

  #============================================================================#
  # Import Data text
  #============================================================================#
  output$import_summary <- renderUI({

    # --- Main title ---
    text_output <- paste0(
      "<h4 style='color:#669999; font-weight:bold; margin-bottom:15px'>",
      "<i class='fa-solid fa-circle-info' style='color:#669999; margin-right:8px;'></i>",
      "Data Import</h4>"
    )

      if (!is.null(rtdataReactive())) {
        polarity <- rtdataReactive()$metadata$polarity
        source <- rtdataReactive()$metadata$source

        # --- CSV ---
        if (input$data_format == "data_file") {
          req(input$path_data_file)
          file_name <- input$path_data_file$name
          file_size <- round(file.info(input$path_data_file$datapath)$size / (1024^2), 2)

          csv_table <- paste0(
            "<hr style='border: 1.5px solid #95a5a6;'>",
            "<p><i class='fa-solid fa-square-check' style='margin-left:20px; color:#95a5a6;'></i><b> File</b></p>",

            "<div style='margin-left:20px; margin-right:80px; border:1px solid #ddd;
            border-radius:6px; padding:5px; background-color:#ffffff; font-size:13px;
            color:#2C3E50; max-width:100%; overflow-x:auto; overflow-y:hidden;'>",

            "<p style='margin-left:20px; margin-top:10px; color:#669999; font-style:italic;'><b> File successfully imported and validated.</b></p>",
            "<hr style='margin-top:15px; margin-bottom:10px; border-top:1px dashed #ccc;'>",

            "<table style='margin-left:50px; border-collapse: collapse; width: auto;'>",

            "<tr style='border-bottom: 1px solid #ccc;'>",
            "<td style='padding:6px 20px 6px 10px; font-weight:bold;'> Name</td>",
            "<td style='padding:6px 10px;'>", file_name, "</td>",
            "</tr>",

            "<tr style='border-bottom: 1px solid #ccc;'>",
            "<td style='padding:6px 20px 6px 10px; font-weight:bold;'> Size</td>",
            "<td style='padding:6px 10px;'>", file_size, " MB</td>",
            "</tr>",

            "<tr style='border-bottom: 1px solid #ccc;'>",
            "<td style='padding:6px 20px 6px 10px; font-weight:bold;'> Source</td>",
            "<td style='padding:6px 10px;'>", source, "</td>",
            "</tr>",

            "<tr>",
            "<td style='padding:6px 20px 6px 10px; font-weight:bold;'> Polarity</td>",
            "<td style='padding:6px 10px;'>", polarity, "</td>",
            "</tr>",

            "</table>",

            "</div>"
          )
          text_output <- paste0(text_output, csv_table)

          # --- TXT ---
        } else if (input$data_format == "data_rdata") {
          req(input$path_data_rdata)
          file_name <- input$path_data_rdata$name
          file_size <- round(file.info(input$path_data_rdata$datapath)$size / (1024^2), 2)

          csv_table <- paste0(
            "<hr style='border: 1.5px solid #95a5a6;'>",
            "<p><i class='fa-solid fa-square-check' style='margin-left:20px; color:#95a5a6;'></i><b> RData File</b></p>",

            "<div style='margin-left:20px; margin-right:80px; border:1px solid #ddd;
            border-radius:6px; padding:5px; background-color:#ffffff; font-size:13px;
            color:#2C3E50; max-width:100%; overflow-x:auto; overflow-y:hidden;'>",

            "<p style='margin-left:20px; margin-top:10px; color:#669999; font-style:italic;'><b> File successfully imported and validated.</b></p>",
            "<hr style='margin-top:15px; margin-bottom:10px; border-top:1px dashed #ccc;'>",

            "<table style='margin-left:50px; border-collapse: collapse; width: auto;'>",

            "<tr style='border-bottom: 1px solid #ccc;'>",
            "<td style='padding:6px 20px 6px 10px; font-weight:bold;'> Name</td>",
            "<td style='padding:6px 10px;'>", file_name, "</td>",
            "</tr>",

            "<tr style='border-bottom: 1px solid #ccc;'>",
            "<td style='padding:6px 20px 6px 10px; font-weight:bold;'> Size</td>",
            "<td style='padding:6px 10px;'>", file_size, " MB</td>",
            "</tr>",

            "<tr style='border-bottom: 1px solid #ccc;'>",
            "<td style='padding:6px 20px 6px 10px; font-weight:bold;'> Source</td>",
            "<td style='padding:6px 10px;'>", source, "</td>",
            "</tr>",

            "<tr>",
            "<td style='padding:6px 20px 6px 10px; font-weight:bold;'> Polarity</td>",
            "<td style='padding:6px 10px;'>", polarity, "</td>",
            "</tr>",

            "</table>",

            "</div>"
          )
          text_output <- paste0(text_output, csv_table)

          # --- TXT ---
        } else if (input$data_format == "data_renv") {
          obj_name <- input$name_data_renv
          obj_size <- round(object.size(rtdataReactive()) / (1024^2), 2)

          renv_table <- paste0(
            "<hr style='border: 1.5px solid #95a5a6;'>",
            "<p><i class='fa-solid fa-square-check' style='margin-left:20px; color:#95a5a6;'></i><b> R environment object</b></p>",

            "<div style='margin-left:20px; margin-right:80px; border:1px solid #ddd;
            border-radius:6px; padding:5px; background-color:#ffffff; font-size:13px;
            color:#2C3E50; max-width:100%; overflow-x:auto; overflow-y:hidden;'>",

            "<p style='margin-left:20px; margin-top:10px; color:#669999; font-style:italic;'><b> Object successfully loaded and validated.</b></p>",
            "<hr style='margin-top:15px; margin-bottom:10px; border-top:1px dashed #ccc;'>",

            "<table style='margin-left:50px; border-collapse: collapse; width: auto;'>",

            "<tr style='border-bottom: 1px solid #ccc;'>",
            "<td style='padding:6px 20px 6px 10px; font-weight:bold;'>️ Name</td>",
            "<td style='padding:6px 10px;'>", obj_name, "</td>",
            "</tr>",

            "<tr style='border-bottom: 1px solid #ccc;'>",
            "<td style='padding:6px 20px 6px 10px; font-weight:bold;'> Size</td>",
            "<td style='padding:6px 10px;'>", obj_size, " MB</td>",
            "</tr>",

            "<tr style='border-bottom: 1px solid #ccc;'>",
            "<td style='padding:6px 20px 6px 10px; font-weight:bold;'> Source</td>",
            "<td style='padding:6px 10px;'>", source, "</td>",
            "</tr>",

            "<tr>",
            "<td style='padding:6px 20px 6px 10px; font-weight:bold;'> Polarity</td>",
            "<td style='padding:6px 10px;'>", polarity, "</td>",
            "</tr>",

            "</table>",

            "</div>"
          )

          text_output <- paste0(text_output, renv_table)
        }

      } else {
        # --- No import yet ---
        if (input$data_format == "data_file") {

          text_output <- paste0(
            text_output,
            "<hr style='border: 1.5px solid #95a5a6;'>",
            "<p><i class='fa-regular fa-square' style='margin-left:20px; color:gray;'></i><b> Peaklist File</b></p>",

            "<div style='margin-left:20px; margin-right:80px; border:1px solid #ddd;
            border-radius:6px; padding:5px; background-color:#ffffff; font-size:13px;
            color:#2C3E50; max-width:100%; overflow-x:auto; overflow-y:hidden;'>",

            "<p style='margin-left:20px; margin-top:10px; font-style:italic; font-size:13px;'>",
            "No data has been imported yet. Please upload a valid <b>.csv</b>, <b>.xlsx</b> or <b>.txt</b> file.</p>",

            "<hr style='margin-left:40px; margin-right:40px; border-top:1px dashed #ccc;'>",
            "<p style='margin-left:20px; margin-top:10px; font-style:italic; font-size:13px;'>",
            "The file may originate from one of the following sources:</p>",

            "<ul style='margin-left:40px; font-style:italic; font-size:13px;'>",
            "<li><b><span style='color:#669999;'>LipidMS</span></b>: Use the <i>features</i> table from one of the following objects:<br>
            &nbsp;&nbsp;&nbsp;&nbsp;- <b><span style='color:#2C3E50;'>msbatch</span></b>: <span style='color:#2C3E50;'>msbatch$features</span><br>
            &nbsp;&nbsp;&nbsp;&nbsp;- <b><span style='color:#2C3E50;'>msobject</span></b>: <span style='color:#2C3E50;'>msobject[[sample]]$annotation$annotatedPeaklist</span>
            </li>",
            "<hr style='margin-left:40px; margin-right:40px; border-top:1px dashed #ccc;'>",
            "<li><b><span style='color:#669999;'>MS-DIAL</span></b>: Export the <i>Alignment Results</i> in <b>.csv</b> or <b>.txt</b> format.</li>",
            "<hr style='margin-left:40px; margin-right:40px; border-top:1px dashed #ccc;'>",
            "<li><b><span style='color:#669999;'>Other</span></b>: Provide a custom <b>.csv</b> or <b>.txt</b> file that includes the following columns:</li>",
            "</ul>",

            "<table class='no-hover' style='margin-left:120px; margin-bottom:10px; border-collapse: collapse; font-size:13px; font-style:italic;'>",
            "<thead>",
            "<tr style='background-color:#ecf0f1;'>",
            "<th style='border: 1px solid #ddd; padding: 6px; white-space: nowrap; text-align:center; text-decoration: underline;'>mz</th>",
            "<th style='border: 1px solid #ddd; padding: 6px; white-space: nowrap; text-align:center; text-decoration: underline;'>RTsec or RTmin</th>",
            "<th style='border: 1px solid #ddd; padding: 6px; white-space: nowrap; text-align:center; text-decoration: underline;'>lipidname</th>",
            "<th style='border: 1px solid #ddd; padding: 6px; white-space: nowrap; text-align:center; text-decoration: underline;'>adduct</th>",
            "<th style='border: 1px solid #ddd; padding: 6px; white-space: nowrap; text-align:center; text-decoration: underline;'>confidenceLevel</th>",
            "</tr>",
            "</thead>",
            "<tbody>",
            "<tr>",
            "<td style='border: 1px solid #ddd; padding: 6px; text-align:center'>Mass-to-charge<br>ratio</td>",
            "<td style='border: 1px solid #ddd; padding: 6px; text-align:center'>Retention time<br>(in seconds or minutes)</td>",
            "<td style='border: 1px solid #ddd; padding: 6px; text-align:center'>Name of the<br>lipid</td>",
            "<td style='border: 1px solid #ddd; padding: 6px; text-align:center'>Ion type<br>(e.g., [M+H]+)</td>",
            "<td style='border: 1px solid #ddd; padding: 6px; text-align:center'>Reliability of the annotation<br>(e.g., FA position)</td>",
            "</tr>",
            "</tbody>",
            "</table>",

            "<p style='margin-left:120px; color:#B03A2E; font-size:12px; font-style:italic;'>",
            "Ensure the file contains no empty rows, that all required column headers are present, ",
            "and that the correct column and adduct delimiters are selected to ensure proper import.",
            "</p>",
            "</div>"
          )

        } else if (input$data_format == "data_rdata") {

          text_output <- paste0(
            text_output,
            "<hr style='border: 1.5px solid #95a5a6;'>",
            "<p><i class='fa-regular fa-square' style='margin-left:20px; color:gray;'></i><b> RData File</b></p>",

            "<div style='margin-left:20px; margin-right:80px; border:1px solid #ddd;
            border-radius:6px; padding:5px; background-color:#ffffff; font-size:13px;
            color:#2C3E50; max-width:100%; overflow-x:auto; overflow-y:hidden;'>",

            "<p style='margin-left:20px; margin-top:10px; font-style:italic; font-size:13px;'>",
            "No RData File has been imported yet. Please upload a valid <b>.RData</b> or <b>.rda.gz</b> file.</p>",

            "<hr style='margin-left:40px; margin-right:40px; border-top:1px dashed #ccc;'>",
            "<p style='margin-left:20px; margin-top:10px; font-style:italic; font-size:13px;'>",
            "The file must originate from <b><span style='color:#669999;'>LipidRT</span></b>.</p>",
            "<ul style='margin-left:40px; font-style:italic; font-size:13px;'>",
            "<li>",
            "Import the full <b><span style='color:#2C3E50'>rtdata</span></b> object directly.<br>",
            "</li>",
            "</ul>",

            "</div>"
          )

        } else if (input$data_format == "data_renv") {

          text_output <- paste0(
            text_output,
            "<hr style='border: 1.5px solid #95a5a6;'>",
            "<p><i class='fa-regular fa-square' style='margin-left:20px; color:gray;'></i><b> R environment object</b></p>",

            "<div style='margin-left:20px; margin-right:80px; border:1px solid #ddd;
            border-radius:6px; padding:5px; background-color:#ffffff; font-size:13px;
            color:#2C3E50; max-width:100%; overflow-x:auto; overflow-y:hidden;'>",

            "<p style='margin-left:20px; margin-top:10px; font-style:italic; font-size:13px;'>",
            "To import data from the R environment, please enter the <b>name of the variable</b> that contains the data table in the input field and click the <b>Import</b> button.</p>",

            "<hr style='margin-left:40px; margin-right:40px; border-top:1px dashed #ccc;'>",
            "<p style='margin-left:20px; margin-top:10px; font-style:italic; font-size:13px;'>",
            "The object may originate from one of the following sources:</p>",

            "<ul style='margin-left:40px; font-style:italic; font-size:13px;'>",
            "<li><b><span style='color:#669999;'>LipidRT</span></b>: Load the <b>rtdata object</b> returned by the <i>readrtdata function</i>. The object may already be filtered and predicted.</li>",
            "<hr style='margin-left:40px; margin-right:40px; border-top:1px dashed #ccc;'>",
            "<li><b><span style='color:#669999;'>LipidMS</span></b>: Use the <i>features</i> table from one of the following objects:<br>
            &nbsp;&nbsp;&nbsp;&nbsp;- <b><span style='color:#2C3E50;'>msbatch</span></b>: <span style='color:#2C3E50;'>msbatch$features</span><br>
            &nbsp;&nbsp;&nbsp;&nbsp;- <b><span style='color:#2C3E50;'>msobject</span></b>: <span style='color:#2C3E50;'>msobject[[sample]]$annotation$annotatedPeaklist</span><br>
            Assign the table to a variable before importing.</li>",
            "<hr style='margin-left:40px; margin-right:40px; border-top:1px dashed #ccc;'>",
            "<li><b><span style='color:#669999;'>MS-DIAL</span></b>: Load the <i>Alignment Results</i> as a data frame from a <b>.csv</b> or <b>.txt</b> file and assign it to a variable.</li>",
            "<hr style='margin-left:40px; margin-right:40px; border-top:1px dashed #ccc;'>",
            "<li><b><span style='color:#669999;'>Other</span></b>: Provide a custom <b>data frame</b> with the following columns:</li>",
            "</ul>",

            "<table class='no-hover' style='margin-left:120px; margin-bottom:10px; border-collapse: collapse; font-size:13px; font-style:italic;'>",
            "<thead>",
            "<tr style='background-color:#ecf0f1;'>",
            "<th style='border: 1px solid #ddd; padding: 6px; white-space: nowrap; text-align:center; text-decoration: underline;'>mz</th>",
            "<th style='border: 1px solid #ddd; padding: 6px; white-space: nowrap; text-align:center; text-decoration: underline;'>RT</th>",
            "<th style='border: 1px solid #ddd; padding: 6px; white-space: nowrap; text-align:center; text-decoration: underline;'>lipidname</th>",
            "<th style='border: 1px solid #ddd; padding: 6px; white-space: nowrap; text-align:center; text-decoration: underline;'>adduct</th>",
            "<th style='border: 1px solid #ddd; padding: 6px; white-space: nowrap; text-align:center; text-decoration: underline;'>confidenceLevel</th>",
            "</tr>",
            "</thead>",
            "<tbody>",
            "<tr>",
            "<td style='border: 1px solid #ddd; padding: 6px; text-align:center;'>Mass-to-charge<br>ratio</td>",
            "<td style='border: 1px solid #ddd; padding: 6px; text-align:center;'>Retention time<br>(in seconds)</td>",
            "<td style='border: 1px solid #ddd; padding: 6px; text-align:center;'>Name of the<br>lipid</td>",
            "<td style='border: 1px solid #ddd; padding: 6px; text-align:center;'>Ion type<br>(e.g., [M+H]+)</td>",
            "<td style='border: 1px solid #ddd; padding: 6px; text-align:center;'>Reliability of the annotation<br>(e.g., FA position)</td>",
            "</tr>",
            "</tbody>",
            "</table>",

            "<p style='margin-left:120px; color:#B03A2E; font-size:12px; font-style:italic;'>",
            "Ensure the file contains no empty rows, that all required column headers are present, ",
            "and that the correct column and adduct delimiters are selected to ensure proper import.",
            "</p>",

            "<hr style='margin-left:40px; margin-right:40px; border-top:1px dashed #ccc;'>",
            "<p style='margin-left:20px; margin-top:10px; font-style:italic; font-size:13px;'>",
            "You must provide a variable name (e.g., <code>mydata</code>) that already exists in the R environment and contains the required table. Expressions like <code>msbatch$features</code> are not valid — instead, assign them to a variable first (e.g., <code>mydata &lt;- msbatch$features</code>).</p>",

            "</b>",
            "</div>"
          )

        }
      }


    # Lipid classes (default or custom) ---
    if (!is.null(input$type_lipidclasses)) {

      if (is.null(rtdataReactive()$metadata$polarity)) {

        if (input$data_format == "data_file") {
          data_format <- paste0(
            "The Peaklist File has not been uploaded yet. Upload a valid <b>.csv</b>, <b>.xlsx</b>, or <b>.txt</b> file.<br>",
            "&nbsp;&nbsp;&nbsp;&nbsp; - <b>Polarity is missing</b>"
          )
        } else if (input$data_format == "data_rdata") {
          data_format <- paste0(
            "The RData File has not been uploaded yet. Upload a valid <b>.RData</b> or <b>.rda.gz</b> file.<br>",
            "&nbsp;&nbsp;&nbsp;&nbsp; - <b>Polarity is missing</b>"
          )
        } else if (input$data_format == "data_renv") {
          data_format <- paste0(
            "The data from the R environment has not been uploaded yet. Import a valid <b>object from the R environment</b><br>",
            "&nbsp;&nbsp;&nbsp;&nbsp; - <b>Polarity is missing</b>"
          )
        }

        text_output <- paste0(
          text_output,
          "<hr style='border: 1.5px solid #95a5a6;'>",
          "<p><i class='fa-regular fa-square' style='margin-left:20px; color:gray;'></i><b> Lipid classes for analysis</b></p>"
        )

      } else {

        if (input$data_format == "data_file") {
          data_format <- paste0(
            "The Peaklist File has been successfully uploaded.<br>",
            "&nbsp;&nbsp;&nbsp;&nbsp; - <b>Polarity detected</b>: ",
            rtdataReactive()$metadata$polarity
          )
        } else if (input$data_format == "data_rdata") {
          data_format <- paste0(
            "The RData File has been successfully uploaded.<br>",
            "&nbsp;&nbsp;&nbsp;&nbsp; - <b>Polarity detected</b>: ",
            rtdataReactive()$metadata$polarity
          )
        } else if (input$data_format == "data_renv") {
          data_format <- paste0(
            "The data from the R environment has been successfully uploaded.<br>",
            "&nbsp;&nbsp;&nbsp;&nbsp; - <b>Polarity detected</b>: ",
            rtdataReactive()$metadata$polarity
          )
        }

        text_output <- paste0(
          text_output,
          "<hr style='border: 1.5px solid #95a5a6;'>",
          "<p><i class='fa-solid fa-square-check' style='margin-left:20px; color:#95a5a6;'></i><b> Lipid classes for analysis</b></p>"
        )
      }


      # --- 1: data_format ---
      text_output_dataformat <- paste0(
        "<div style='margin-left:20px; margin-right:80px; border:1px solid #ddd;",
        "border-radius:6px; padding:5px; background-color:#ffffff; font-size:13px;",
        "color:#2C3E50; max-width:100%; overflow-x:auto; overflow-y:hidden;'>",

        "<p style='margin-left:20px; margin-top:10px; font-style:italic; font-size:13px;'>",
        data_format,
        ".</p>",

        "<hr style='border-top:1px dashed #ccc;'>"
      )

      text_output <- paste0(
        text_output,
        text_output_dataformat
      )


      # --- 2: adducts Table ---
      if (input$type_adductsTable == "adductsTable_default") {
        text_type_adductsTable <- "default"
      } else if (input$type_adductsTable == "adductsTable_custom") {
        text_type_adductsTable <- "custom"
      } else {
        text_type_adductsTable <- "unknown"
      }

      # unique ID
      modal_id <- paste0("adductsModal_", as.integer(runif(1, 1000, 9999)))

      # adductsTable table
      adductsTable <- adductsTableReactive()
      if (is.null(adductsTable) || nrow(adductsTable) == 0) {
        adducts_text <- "<div style='margin-left:30px; color:gray; font-style:italic;'>No adductsTable selected.</div>"
      } else {
        adducts_text <- "<table class='no-hover' style='border-collapse: collapse; margin:15px auto; font-size:13px; color:#2C3E50; width:auto; max-width:100%; table-layout:auto; box-shadow: 0 2px 5px rgba(0,0,0,0.05); border-radius:4px;'>"

        # Header
        adducts_text <- paste0(
          adducts_text,
          "<thead><tr style='background-color:#ecf0f1;'>",
          "<th style='border:1px solid #ddd; padding:8px 12px; text-align:center; font-weight:bold; text-decoration:underline;'>adduct</th>",
          "<th style='border:1px solid #ddd; padding:8px 12px; text-align:center; font-weight:bold; text-decoration:underline;'>mdiff</th>",
          "<th style='border:1px solid #ddd; padding:8px 12px; text-align:center; font-weight:bold; text-decoration:underline;'>charge</th>",
          "<th style='border:1px solid #ddd; padding:8px 12px; text-align:center; font-weight:bold; text-decoration:underline;'>n</th>",
          "</tr></thead><tbody>"
        )

        # Rows
        for (i in seq_len(nrow(adductsTable))) {
          adducts_text <- paste0(
            adducts_text,
            "<tr style='transition:background-color 0.2s; cursor:default;' onmouseover='this.style.background=\"#f8f8f8\";' onmouseout='this.style.background=\"white\";'>",
            "<td style='border:1px solid #ddd; text-align:center; max-width:350px; white-space:nowrap; overflow:hidden; text-overflow:ellipsis; padding:6px 12px;' title='", adductsTable$adduct[i], "'><b>", adductsTable$adduct[i], "</b></td>",
            "<td style='border:1px solid #ddd; text-align:center; max-width:350px; white-space:nowrap; overflow:hidden; text-overflow:ellipsis; padding:6px 12px;' title='", adductsTable$mdiff[i], "'>", adductsTable$mdiff[i], "</td>",
            "<td style='border:1px solid #ddd; text-align:center; max-width:350px; white-space:nowrap; overflow:hidden; text-overflow:ellipsis; padding:6px 12px;' title='", adductsTable$charge[i], "'>", adductsTable$charge[i], "</td>",
            "<td style='border:1px solid #ddd; text-align:center; max-width:350px; white-space:nowrap; overflow:hidden; text-overflow:ellipsis; padding:6px 12px;' title='", adductsTable$n[i], "'>", adductsTable$n[i], "</td>",
            "</tr>"
          )
        }

        adducts_text <- paste0(adducts_text, "</tbody></table>")
      }


      text_output_adducts <- paste0(
        "<div style='display:flex; align-items:center; background-color:#ecf0f1; padding:4px 8px; border-radius:4px; margin-left:20px; font-style:italic; font-size:13px;'>",
        "<p style='margin:0; font-size:14px; color:black;'>",
        "<b style='text-decoration: underline;'>Adducts Table:</b> ",
        "<span style='color:black;'>(", text_type_adductsTable, " table)</span>",
        "</p>",
        "<button onclick='document.getElementById(\"", modal_id, "\").style.display=\"flex\"' ",
        "style='margin-left:auto; margin-right:80px; padding:4px 10px; font-size:13px; cursor:pointer; border-radius:3px; background-color:#f9f9f9; color:black; border:1px solid #ccc;'>",
        "👁 Show current Adducts Table",
        "</button>",
        "</div>",

        # Text
        "<p style='margin-top:10px; margin-left:90px; color:#2C3E50; font-style:italic; font-size:14px;'>",
        "You can modify the current adducts table by importing a <b>csv file</b> or a <b>data frame from the global environment</b>.<br>",
        "The Adducts Table must follow the following format:",
        "</p>",

        # Table
        "<table class='no-hover' style='margin-left:120px; margin-bottom:10px; border-collapse: collapse; font-size:13px; font-style:italic;'>",
        "<tr style='background-color:#f9f9f9;'>",
        "<th style='border:1px solid #ddd; padding:6px; text-align:center; text-decoration:underline;'>adduct</th>",
        "<th style='border:1px solid #ddd; padding:6px; text-align:center; text-decoration:underline;'>mdiff</th>",
        "<th style='border:1px solid #ddd; padding:6px; text-align:center; text-decoration:underline;'>charge</th>",
        "<th style='border:1px solid #ddd; padding:6px; text-align:center; text-decoration:underline;'>n</th>",
        "</tr>",
        "<tbody>",
        "<tr>",
        "<td style='border:1px solid #ddd; padding:6px; text-align:center;'>Adduct notation<br>(e.g. M+H, M+Na, M-H)</td>",
        "<td style='border:1px solid #ddd; padding:6px; text-align:center;'>Mass difference added to the molecule<br>(e.g. 1.0073 for H)</td>",
        "<td style='border:1px solid #ddd; padding:6px; text-align:center;'>Ion charge state<br>(e.g. 1, -1)</td>",
        "<td style='border:1px solid #ddd; padding:6px; text-align:center;'>Number of molecules in the adduct<br>(e.g. 1 for M+H, 2 for 2M+H)</td>",
        "</tr>",
        "</tbody>",
        "</table>",
        "<p style='margin-left:120px; color:#B03A2E; font-size:12px; font-style:italic;'>
          Ensure the file contains no empty rows and that all required column headers are present.
        The correct column delimiters must be used for proper import.
        </p>",

        # Modal
        "<div id='", modal_id, "' style='display:none; position:fixed; z-index:9999; left:0; top:0; width:100%; height:100%; background-color: rgba(0,0,0,0.5); justify-content:center; align-items:center;'>",
        "<div style='background-color:white; padding:20px; border-radius:6px; min-width:50%; max-width:80%; max-height:80%; overflow:auto; position:relative; display:inline-block;'>",
        "<span onclick='document.getElementById(\"", modal_id, "\").style.display=\"none\"' ",
        "style='position:absolute; top:10px; right:20px; color:#333; font-size:24px; font-weight:bold; cursor:pointer;'>&times;</span>",
        adducts_text,
        "</div>",

        # JS JS to close the modal when clicking outside it
        "<script>
  var modal_adducts = document.getElementById('", modal_id, "');
  window.addEventListener('click', function(event) {
    if(event.target == modal_adducts) {
      modal_adducts.style.display = 'none';
    }
  });
</script>",

        "</div>"
      )

      text_output <- paste0(
        text_output,
        text_output_adducts,
        "<hr style='border-top:1px dashed #ccc;'>"
        )

      # --- 3: filedbs Table ---
      if (input$type_filedbs == "filedbs_default") {
        text_type_filedbs <- "default"
      } else if (input$type_filedbs == "filedbs_custom") {
        text_type_filedbs <- "custom"
      } else {
        text_type_filedbs <- "unknown"
      }

      # unique ID
      modal_id_filedbs <- paste0("filedbsModal_", as.integer(runif(1, 1000, 9999)))

      filedbs <- filedbsReactive()
      if (is.null(filedbs) || nrow(filedbs) == 0) {
        filedbs_text <- "<div style='margin-left:30px; color:gray; font-style:italic;'>No File dbs Table selected.</div>"
      } else {
        filedbs_text <- "<table class='no-hover' style='border-collapse: collapse; margin:15px auto; font-size:13px; color:#2C3E50; width:auto; max-width:100%; table-layout:auto; box-shadow: 0 2px 5px rgba(0,0,0,0.05); border-radius:4px;'>"

        # Header
        filedbs_text <- paste0(
          filedbs_text,
          "<thead><tr style='background-color:#ecf0f1;'>",
          "<th style='border:1px solid #ddd; padding:8px 12px; text-align:center; font-weight:bold; text-decoration:underline;'>classes</th>",
          "<th style='border:1px solid #ddd; padding:8px 12px; text-align:center; font-weight:bold; text-decoration:underline;'>adducts_pos</th>",
          "<th style='border:1px solid #ddd; padding:8px 12px; text-align:center; font-weight:bold; text-decoration:underline;'>adducts_neg</th>",
          "<th style='border:1px solid #ddd; padding:8px 12px; text-align:center; font-weight:bold; text-decoration:underline;'>fas</th>",
          "<th style='border:1px solid #ddd; padding:8px 12px; text-align:center; font-weight:bold; text-decoration:underline;'>sph</th>",
          "<th style='border:1px solid #ddd; padding:8px 12px; text-align:center; font-weight:bold; text-decoration:underline;'>alkyl</th>",
          "<th style='border:1px solid #ddd; padding:8px 12px; text-align:center; font-weight:bold; text-decoration:underline;'>vinyl</th>",
          "</tr></thead><tbody>"
        )

        # Rows
        lipidclasses <- unique(filedbs$classes)
        for (lipid in lipidclasses) {
          adductspos_text <- gsub("\\s*([|;,])\\s*", " \\1 ", filedbs[["adducts_pos"]][filedbs$classes == lipid])
          adductsneg_text <- gsub("\\s*([|;,])\\s*", " \\1 ", filedbs[["adducts_neg"]][filedbs$classes == lipid])
          fas_text <- gsub("\\s*([|;,])\\s*", " \\1 ", filedbs$fas[filedbs$classes == lipid])
          sph_text <- gsub("\\s*([|;,])\\s*", " \\1 ", filedbs$sph[filedbs$classes == lipid])
          alkyl_text <- gsub("\\s*([|;,])\\s*", " \\1 ", filedbs$alkyl[filedbs$classes == lipid])
          vinyl_text <- gsub("\\s*([|;,])\\s*", " \\1 ", filedbs$vinyl[filedbs$classes == lipid])

          filedbs_text <- paste0(
            filedbs_text,
            "<tr style='transition:background-color 0.2s; cursor:default;' onmouseover='this.style.background=\"#f8f8f8\";' onmouseout='this.style.background=\"white\";'>",
            "<td style='border:1px solid #ddd; text-align:center; max-width:350px; white-space:normal; overflow-wrap: break-word; word-break:break-word; padding:6px 12px;' title='", lipid, "'><b>", lipid, "</b></td>",
            "<td style='border:1px solid #ddd; text-align:center; max-width:350px; white-space:normal; overflow-wrap: break-word; word-break:break-word; padding:6px 12px;' title='", adductspos_text, "'>", adductspos_text, "</td>",
            "<td style='border:1px solid #ddd; text-align:center; max-width:350px; white-space:normal; overflow-wrap: break-word; word-break:break-word; padding:6px 12px;' title='", adductsneg_text, "'>", adductsneg_text, "</td>",
            "<td style='border:1px solid #ddd; text-align:center; max-width:350px; white-space:normal; overflow-wrap: break-word; word-break:break-word; padding:6px 12px;' title='", fas_text, "'>", fas_text, "</td>",
            "<td style='border:1px solid #ddd; text-align:center; max-width:350px; white-space:normal; overflow-wrap: break-word; word-break:break-word; padding:6px 12px;' title='", sph_text, "'>", sph_text, "</td>",
            "<td style='border:1px solid #ddd; text-align:center; max-width:350px; white-space:normal; overflow-wrap: break-word; word-break:break-word; padding:6px 12px;' title='", alkyl_text, "'>", alkyl_text, "</td>",
            "<td style='border:1px solid #ddd; text-align:center; max-width:350px; white-space:normal; overflow-wrap: break-word; word-break:break-word; padding:6px 12px;' title='", vinyl_text, "'>", vinyl_text, "</td>",
            "</tr>"
          )
        }
        filedbs_text <- paste0(filedbs_text, "</tbody></table>")
      }

      text_output_filedbs <- paste0(
        "<div style='display:flex; align-items:center; background-color:#ecf0f1; padding:4px 8px; border-radius:4px; margin-left:20px; font-style:italic; font-size:13px;'>",
        "<p style='margin:0; font-size:14px; color:black;'>",
        "<b style='text-decoration: underline;'>File dbs Table:</b> ",
        "<span style='color:black;'>(", text_type_filedbs, " table)</span>",
        "</p>",
        "<button onclick='document.getElementById(\"", modal_id_filedbs, "\").style.display=\"flex\"' ",
        "style='margin-left:auto; margin-right:80px; padding:4px 10px; font-size:13px; cursor:pointer; border-radius:3px; background-color:#f9f9f9; color:black; border:1px solid #ccc;'>",
        "👁 Show current File dbs Table",
        "</button>",
        "</div>",

        # Text
        "<p style='margin-top:10px; margin-left:90px; color:#2C3E50; font-style:italic; font-size:14px;'>",
        "You can modify the current file dbs table by importing a <b>csv file</b> or a <b>data frame from the global environment</b>.<br>",
        "The file dbs Table must follow the following format:",
        "</p>",

        # Table
        "<table class='no-hover' style='margin-left:120px; margin-bottom:10px; border-collapse: collapse; font-size:13px; font-style:italic;'>",
        "<tr style='background-color:#f9f9f9;'>",
        "<th style='border:1px solid #ddd; padding:6px; text-align:center; text-decoration:underline;'>classes</th>",
        "<th style='border:1px solid #ddd; padding:6px; text-align:center; text-decoration:underline;'>adducts_pos</th>",
        "<th style='border:1px solid #ddd; padding:6px; text-align:center; text-decoration:underline;'>adducts_neg</th>",
        "<th style='border:1px solid #ddd; padding:6px; text-align:center; text-decoration:underline;'>fas</th>",
        "<th style='border:1px solid #ddd; padding:6px; text-align:center; text-decoration:underline;'>sph</th>",
        "<th style='border:1px solid #ddd; padding:6px; text-align:center; text-decoration:underline;'>alkyl</th>",
        "<th style='border:1px solid #ddd; padding:6px; text-align:center; text-decoration:underline;'>vinyl</th>",
        "</tr>",
        "<tbody>",
        "<tr>",
        "<td style='border:1px solid #ddd; padding:6px; text-align:center;'>Lipid class<br>(e.g. PC)</td>",
        "<td style='border:1px solid #ddd; padding:6px; text-align:center;'>Positive-mode adducts<br>(e.g. M+H|M+Na)</td>",
        "<td style='border:1px solid #ddd; padding:6px; text-align:center;'>Negative-mode adducts<br>(e.g. M-H|M+CH3COO)</td>",
        "<td style='border:1px solid #ddd; padding:6px; text-align:center;'>Fatty acid chains<br>(e.g. 8:0|10:0|12:0|...)</td>",
        "<td style='border:1px solid #ddd; padding:6px; text-align:center;'>Sphingoid bases<br>(e.g. 16:0|18:1|...)</td>",
        "<td style='border:1px solid #ddd; padding:6px; text-align:center;'>Alkyl chains<br>(e.g. 16:0|18:1|...)</td>",
        "<td style='border:1px solid #ddd; padding:6px; text-align:center;'>Vinyl chains<br>(e.g. 16:0|18:0)</td>",
        "</tr>",
        "</tbody>",
        "</table>",
        "<p style='margin-left:120px; color:#B03A2E; font-size:12px; font-style:italic;'>
          Ensure the file contains no empty rows and that all required column headers are present. The correct column delimiters must be used for proper import.<br>
          The columns <b>adducts_pos</b> and <b>adducts_neg</b> can both be provided, or only the one corresponding to the analysis polarity.
          However, the column matching the analysis polarity is mandatory.<br>
          <b>Note:</b> Only the adducts listed in the Adducts Table will be included.
          If you need to add extra adducts, modify the Adducts Table accordingly to ensure they are considered in the analysis.
        </p>",

        # Modal - window
        "<div id='", modal_id_filedbs, "' style='display:none; position:fixed; z-index:9999; left:0; top:0; width:100%; height:100%; background-color: rgba(0,0,0,0.5); justify-content:center; align-items:center;'>",
        "<div style='background-color:white; padding:20px; border-radius:6px; min-width:50%; max-width:80%; max-height:80%; overflow:auto; position:relative; display:inline-block;'>",
        "<span onclick='document.getElementById(\"", modal_id_filedbs, "\").style.display=\"none\"' ",
        "style='position:absolute; top:10px; right:20px; color:#333; font-size:24px; font-weight:bold; cursor:pointer;'>&times;</span>",
        filedbs_text,
        "</div>",

        # JS to close the modal when clicking outside it
        "<script>
  var modal_filedbs = document.getElementById('", modal_id_filedbs, "');
  window.addEventListener('click', function(event) {
    if(event.target == modal_filedbs) {
      modal_filedbs.style.display = 'none';
    }
  });
</script>",

        "</div>"
      )

      text_output <- paste0(text_output,
                            text_output_filedbs)



    }

    text_output <- paste0(text_output,
                          "</div>",
                          "<hr style='border: 1.5px solid #95a5a6; margin-top: 20px; margin-bottom: 20px;'>")

    HTML(text_output)
  })


  # output$dmz_predict_ui <- renderUI({
  #
  #   if (isTRUE(input$ppm_predict)) {
  #
  #     sliderInput(
  #       "dmz_predict",
  #       label = NULL,
  #       min = 0,
  #       max = 50,
  #       value = 5,
  #       step = 1
  #     )
  #
  #   } else {
  #
  #     numericInput(
  #       "dmz_predict",
  #       label = NULL,
  #       value = 0.005,
  #       min = 0,
  #       step = 0.001
  #     )
  #
  #   }
  #
  # })
  observeEvent(input$ppm_predict, {

    if (isTRUE(input$ppm_predict)) {

      updateNumericInput(
        session,
        "dmz_predict",
        value = 5,
        min = 0,
        step = 1
      )

    } else {

      updateNumericInput(
        session,
        "dmz_predict",
        value = 0.005,
        min = 0,
        step = 0.001
      )

    }

  }, ignoreInit = TRUE)



  observeEvent(
    list(rtdataReactive(), input$type_lipidclasses,
         input$type_filedbs, input$type_adductsTable), {

    # If default lipidclasses...
    if (input$type_lipidclasses == "lipidclasses_default") {

      filedbspathReactive(NULL)
      adductsTablepathReactive(NULL)

      output$fileadductsTable_upload_ui <- renderUI({
        fileInput("fileadductsTable", NULL,
                  multiple = FALSE,
                  accept = c("text/plain", ".csv"),
                  width = "100%",
                  placeholder = "No file selected")
      })

      output$filedbs_upload_ui <- renderUI({
        fileInput("filedbs", NULL,
                  multiple = FALSE,
                  accept = c("text/plain", ".csv"),
                  width = "100%",
                  placeholder = "No file selected")
      })


      generateadductsmzdbs <- LipidRT::createadductsmzdbs(filedbs = NULL,
                                                          sep = input$sep_filedbs,
                                                          adducts = NULL)

      filedbs <- generateadductsmzdbs$filedbs
      dbs <- generateadductsmzdbs$adductsmzdbs

      adductsTableReactive(LipidRT::adductsTable)
      filedbsReactive(filedbs)
      dbsReactive(dbs)


      select_lipidclasses(analysis_step = "lipidclasses_predict",
                          type_lipidclasses = "lipidclasses_default")


      # If custom filedbs ...
    } else if (input$type_lipidclasses == "lipidclasses_custom") {

      if (input$type_filedbs == "filedbs_default") {
        filedbspathReactive(NULL)

        output$filedbs_upload_ui <- renderUI({
          fileInput("filedbs", NULL,
                    multiple = FALSE,
                    accept = c("text/plain", ".csv"),
                    width = "100%",
                    placeholder = "No file selected")
        })


        generateadductsmzdbs <- LipidRT::createadductsmzdbs(filedbs = NULL,
                                                            sep = input$sep_filedbs,
                                                            adducts = adductsTableReactive())

        filedbs <- generateadductsmzdbs$filedbs
        dbs <- generateadductsmzdbs$adductsmzdbs

        filedbsReactive(filedbs)
        dbsReactive(dbs)


        select_lipidclasses(analysis_step = "lipidclasses_predict",
                            type_lipidclasses = "lipidclasses_default")

        } else if (input$type_filedbs == "filedbs_custom") {

          select_lipidclasses(analysis_step = "lipidclasses_predict",
                              type_lipidclasses = "lipidclasses_custom")
        }


      if (input$type_adductsTable == "adductsTable_default") {
        adductsTableReactive(LipidRT::adductsTable)
        adductsTablepathReactive(NULL)


          output$adductsTable_upload_ui <- renderUI({
            fileInput("fileadductsTable", NULL,
                      multiple = FALSE,
                      accept = c("text/plain", ".csv"),
                      width = "100%",
                      placeholder = "No file selected")
          })

          generateadductsmzdbs <- LipidRT::createadductsmzdbs(filedbs = filedbspathReactive(),
                                                              sep = input$sep_filedbs,
                                                              adducts = adductsTableReactive())

          filedbs <- generateadductsmzdbs$filedbs
          dbs <- generateadductsmzdbs$adductsmzdbs

          filedbsReactive(filedbs)
          dbsReactive(dbs)

      } else if (input$type_adductsTable == "adductsTable_custom") {

        select_lipidclasses(analysis_step = "lipidclasses_predict",
                            type_lipidclasses = "lipidclasses_custom")
      }


    }

      if (!is.null(rtdataReactive())) {
        rtdata <- rtdataReactive()
        rtdata$metadata$prediction$adductsTable <- adductsTableReactive()
        rtdata$metadata$prediction$filedbs <- filedbsReactive()
        rtdata$metadata$prediction$dbs <- dbsReactive()
        rtdataReactive(rtdata)
      }

  })


  observeEvent(input$lipid_class_selector, {
    selectedclass(input$lipid_class_selector)
  })

  output$lipid_class_selector_ui <- renderUI({
    req(dbChoicesReactive())

    db_choices <- dbChoicesReactive()

    lipid_classes <- names(db_choices$filtering)
    lipid_classes <- lipid_classes[!is.na(lipid_classes) & lipid_classes != ""]

    validate(
      need(length(lipid_classes) > 0, "No classes available")
    )

    selected_value <- selectedclass()

    if (is.null(selected_value) || !selected_value %in% lipid_classes) {
      selected_value <- lipid_classes[1]
    }

    selectInput(
      "lipid_class_selector",
      "Select Lipid Class:",
      choices = lipid_classes,
      selected = selected_value,
      width = "70%"
    )
  })



  # Selector + botones
  output$lipid_classes_combined_selector_ui <- renderUI({
    req(dbChoicesReactive())

    db_choices <- dbChoicesReactive()

    lipid_classes <- names(db_choices$filtering)
    lipid_classes <- lipid_classes[!is.na(lipid_classes) & lipid_classes != ""]

    validate(
      need(length(lipid_classes) > 0, "No classes available")
    )

    div(
      style = "
    border: 1px solid #ddd;
    border-radius: 10px;
    padding: 15px;
    margin-top: 10px;
    margin-bottom: 10px;
    background-color: #fafafa;
    box-shadow: 0 1px 3px rgba(0,0,0,0.08);
  ",

    fluidRow(

      # Selector + botones
      column(
        width = 7,

        div(
          style = " display:flex; gap:8px; align-items:flex-start;",

          # selector
          div(
            style = "flex-grow:1;",
            selectInput(
              "lipid_class_selector_combined",
              "Select Lipid Class:",
              choices = lipid_classes,
              selected = lipid_classes[1],
              width = "100%"
            )
          ),

          # botón
          div(
            style = "margin-top:25px;",

            actionButton(
              "add_lipid_class_combined",
              label = "Add",
              icon = icon("plus"),
              style = "height:38px; min-height:38px; padding: 0px 12px; font-weight:bold; font-size:12px;
              background-color:#95a5a6; color:white;  border:none; border-radius:4px;"
              )
          )
        ),

        tags$small(
          style = "color: #777;",
          "Maximum 7 classes can be selected for combined plotting."
        )
      ),


      # Clases seleccionadas
      column(
        width = 5,
        tags$label("Selected Subclasses:"),
        uiOutput("selected_lipid_classes_combined_ui")
      )
    )
    )
  })


  # Añadir clase
  observeEvent(input$add_lipid_class_combined, {

    current_classes <- lipid_classes_combined()
    new_class <- input$lipid_class_selector_combined

    if (!(new_class %in% current_classes) &&
      length(current_classes) < 7) {
      lipid_classes_combined(
        c(current_classes, new_class)
      )
    }
  })


  # Eliminar clase
  # observeEvent(input$remove_lipid_class_combined, {
  #
  #   current_classes <- lipid_classes_combined()
  #   class_to_remove <- input$lipid_class_selector_combined
  #
  #   if (class_to_remove %in% current_classes) {
  #     lipid_classes_combined(
  #       current_classes[current_classes != class_to_remove]
  #     )
  #   }
  # })


  # Mostrar clases seleccionadas de forma bonita
  output$selected_lipid_classes_combined_ui <- renderUI({

    current_classes <- lipid_classes_combined()

    if (length(current_classes) == 0) {
      return(
        tags$span(
          style = "color:#999;",
          "None selected"
        )
      )
    }

    tagList(

      lapply(seq_along(current_classes), function(i) {

        cl <- current_classes[i]

        tags$div(

          style = paste(
            "position:relative;",
            "display:inline-block;",
            "background-color:#ecf0f1;",
            "border:1px solid #2c3e50;",
            "border-radius:14px;",
            "padding:3px 12px;",   # más bajo
            "margin:5px;",
            "font-size:12px;",
            "line-height:1.2;"
          ),

          # texto
          cl,

          # botón remove
          actionLink(
            inputId = paste0("remove_class_combined_", i),
            label = NULL,
            icon = icon("times"),
            style = paste(
              "position:absolute;",
              "top:-7px;",
              "right:-7px;" ,      # sobresale de la caja
              "width:16px;",
              "height:16px;",
              "border-radius:50%;",
              "background:#e74c3c;",
              "border:1px solid white;",
              "display:flex;",
              "align-items:center;",
              "justify-content:center;",
              "font-size:8px;",
              "color:white;",
              "text-decoration:none;",
              "box-shadow:0 1px 3px rgba(0,0,0,0.25);",
              "cursor:pointer;"
            )
          )
        )
      })
    )
  })

  observe({

    current_classes <- lipid_classes_combined()

    lapply(seq_along(current_classes), function(i) {

      observeEvent(
        input[[paste0("remove_class_combined_", i)]],
        {

          class_to_remove <- current_classes[i]

          # quitar de todas
          lipid_classes_combined(
            setdiff(
              lipid_classes_combined(),
              class_to_remove
            )
          )

          selected_lipid_classes_combined(
            setdiff(
              selected_lipid_classes_combined(),
              class_to_remove
            )
          )

        },
        ignoreInit = TRUE
      )

    })

  })

  #============================================================================#
  # Reset rtdata when switching between File, RData, and R environment.
  #============================================================================#
  observeEvent(input$data_format,{
    rtdataReactive(NULL)

    if (input$data_format == "data_file") {
      if (!is.null(input$path_data_rdata)) {
        output$rdata_upload_ui <- renderUI({
          fileInput("path_data_rdata", NULL,
                    multiple = FALSE,
                    accept = c("text/plain", ".RData", ".gz"),
                    width = "100%",
                    placeholder = "No file selected")
        })
      }
      if (input$name_data_renv != "rtdata") {
        updateTextInput(session, "name_data_renv", value = "rtdata")
      }

    } else if (input$data_format == "data_rdata") {
      if (!is.null(input$path_data_file)) {
        output$file_upload_ui <- renderUI({
          fileInput("path_data_file", NULL,
                    multiple = FALSE,
                    accept = c("text/csv", "text/plain", ".csv", ".txt", ".xlsx"),
                    width = "100%",
                    placeholder = "No file selected")
        })
      }
      if (input$name_data_renv != "rtdata") {
        updateTextInput(session, "name_data_renv", value = "rtdata")
      }

    } else if (input$data_format == "data_renv") {
      if (!is.null(input$path_data_file)) {
        output$file_upload_ui <- renderUI({
          fileInput("path_data_file", NULL,
                    multiple = FALSE,
                    accept = c("text/csv", "text/plain", ".csv", ".txt", ".xlsx"),
                    width = "100%",
                    placeholder = "No file selected")
        })
      }
      if (!is.null(input$path_data_rdata)) {
        output$rdata_upload_ui <- renderUI({
          fileInput("path_data_rdata", NULL,
                    multiple = FALSE,
                    accept = c("text/plain", ".RData", ".gz"),
                    width = "100%",
                    placeholder = "No file selected")
        })
      }
    }
  })

  #============================================================================#
  # Import external CSV/TXT data.
  #============================================================================#
  observeEvent(input$path_data_file, {

    # Detect whether it is CSV or TXT based on the selected option
    if (input$data_format == "data_file" && is.null(input$path_data_file)) {
      return(NULL)
    }

    # --- Importing notification ---
    id_loading <- showNotification(
      HTML(paste0("Importing data from a file... <span class='loader'></span>")),
      type = "default",
      duration = NULL
    )

    tryCatch({
      # Read the input data and determine its source: LipidMS, MS-DIAL, or other.
      rtdata <- LipidRT::readrtdata(input$path_data_file$datapath,
                                    sep = input$csv_col_sep,
                                    sheet = input$sheet_file)
      rtdataReactive(rtdata)

      # Update UI
      update_features_ui(session)

      if (!is.null(rtdata$filtering$data)) {
        rtdata$filtering$data <-
          create_other_isobars(rtdata$filtering$data,
                               typeanalysis = "filtering")

        update_filter_ui(session)
        update_filter_combined_ui(session)

      }
      if (!is.null(rtdata$prediction$candidates)) {

        rtdata$prediction$candidates <-
          create_other_isobars(rtdata$prediction$candidates,
                               typeanalysis = "prediction")

        update_prediction_ui(session)
      }

      # Success notification
      showNotification(
        paste0("The data was imported successfully."),
        type = "message",
        duration = NULL
      )

    }, error = function(e) {
      showNotification(
        paste("Error importing the file:", e$message),
        type = "error"
      )
    }, finally = {
      removeNotification(id_loading)
    })
  })

  observeEvent(list(input$csv_col_sep, input$sheet_file), {
    # Only run if a file is already loaded
    if (is.null(rtdataReactive()) & !is.null(input$path_data_file)) {

      # Reload notification
      id_loading <- showNotification(
        HTML(paste0("Re-importing data with updated separators... <span class='loader'></span>")),
        type = "default",
        duration = NULL
      )

      tryCatch({
        # Re-read the file using the new delimiters
        rtdata <- LipidRT::readrtdata(input$path_data_file$datapath,
                                     sep = input$csv_col_sep,
                                     sheet = input$sheet_file)
        rtdataReactive(rtdata)

        # Update UI
        update_features_ui(session)

        if (!is.null(rtdata$filtering$data)) {
          rtdata$filtering$data <-
            create_other_isobars(rtdata$filtering$data,
                                 typeanalysis = "filtering")

          update_filter_ui(session)
          update_filter_combined_ui(session)
        }
        if (!is.null(rtdata$prediction$candidates)) {

          rtdata$prediction$candidates <-
            create_other_isobars(rtdata$prediction$candidates,
                                 typeanalysis = "prediction")

          update_prediction_ui(session)
        }

        # Success notification
        showNotification(
          paste0("The data was re-imported successfully."),
          type = "message",
          duration = NULL
        )

      }, error = function(e) {
        showNotification(
          paste("Error re-importing the file:", e$message),
          type = "error"
        )
      }, finally = {
        removeNotification(id_loading)
      })
    }
  })



  observeEvent(input$path_data_rdata, {

    if (input$data_format == "data_rdata" && is.null(input$path_data_rdata)) {
      return(NULL)
    }

    id_loading <- showNotification(
      HTML("Importing RData object.<br>It may take a few minutes...<span class='loader'></span>"),
      type = "default",
      duration = NULL
    )

    temp_env <- new.env()

    tryCatch({
      load(input$path_data_rdata$datapath, envir = temp_env)

      rtdata_name <- ls(temp_env)[1]
      rtdata <- temp_env[[rtdata_name]]

      rtdata_structure <- LipidRT::checkdatastructure(rtdata)

      if (rtdata_structure$structure == "rtdata") {
        rtdataReactive(rtdata)

      } else if (is.character(rtdata) || is.data.frame(rtdata)) {
        rtdata <- LipidRT::readrtdata(rtdata)
        rtdataReactive(rtdata)

      } else {
        # Show permanent error notification
        showNotification(
          paste0("The object '", rtdata_name, "' does not have the required format."),
          type = "error",
          duration = NULL
        )

        # Reset the fileInput
        output$rdata_upload_ui <- renderUI({
          fileInput("path_data_rdata", NULL,
                    multiple = FALSE,
                    accept = c("text/plain", ".RData", ".gz"),
                    width = "100%",
                    placeholder = "No file selected")
        })

        removeNotification(id_loading)
        return()
      }

      update_features_ui(session)
      if (!is.null(rtdata$filtering$data)) {
        rtdata$filtering$data <-
          create_other_isobars(rtdata$filtering$data,
                               typeanalysis = "filtering")

        update_filter_ui(session)
        update_filter_combined_ui(session)
      }
      if (!is.null(rtdata$prediction$candidates)) {

        rtdata$prediction$candidates <-
          create_other_isobars(rtdata$prediction$candidates,
                               typeanalysis = "prediction")

        update_prediction_ui(session)
      }

      showNotification(
        paste0("The object '", rtdata_name, "' was imported successfully."),
        type = "message",
        duration = NULL
      )

      removeNotification(id_loading)

    }, error = function(e) {
      showNotification(paste("Error importing the object:", e$message), type = "error", duration = NULL)
      output$rdata_upload_ui <- renderUI({
        fileInput("path_data_rdata", NULL,
                  multiple = FALSE,
                  accept = c("text/plain", ".RData", ".gz"),
                  width = "100%",
                  placeholder = "No file selected")
      })
      removeNotification(id_loading)
    })
  })



  #============================================================================#
  # Import rtdata from the global environment into the app.
  #============================================================================#
  observeEvent(input$import_data_renv, {

    rtdata_name <- input$name_data_renv  # Get the name of the loaded object.

    id_loading <- showNotification(
      HTML(paste0("Importing the object ", rtdata_name, " from the global environment... <span class='loader'></span>")),
      type = "default",
      duration = NULL
    )

    # Check if the object exists in the global environment.
    if (!exists(rtdata_name, envir = .GlobalEnv)) {
      showNotification(
        paste0("The object '", rtdata_name, "' does not exist in the global environment."),
        type = "error",
        duration = NULL
      )
      removeNotification(id_loading)
      return()
    }

    # Get the object from the global environment.
    rtdata <- get(rtdata_name, envir = .GlobalEnv)

    # Check if it has the correct structure.
    rtdata_structure <- LipidRT::checkdatastructure(rtdata)

    # If it has the correct structure...
    if (rtdata_structure$structure == "rtdata") {
      rtdataReactive(rtdata)

    } else if (is.character(rtdata) | is.data.frame(rtdata)) {

      # Read the input data and determine its source: LipidMS, MS-DIAL, or other.
      rtdata <- LipidRT::readrtdata(rtdata)
      rtdataReactive(rtdata)

    } else {
      # If it doesn’t have the correct structure, it notifies you and stops.
      showNotification(
        paste0("The object '", rtdata_name, "' does not have the required format."),
        type = "error"
      )
      removeNotification(id_loading)
      return()
    }

    # Update the interface with the results.
    update_features_ui(session)

    # Update the Filter and Prediction Results interface if they exist.
    if (!is.null(rtdata$filtering$data)) { # Filter
      rtdata$filtering$data <-
        create_other_isobars(rtdata$filtering$data,
                             typeanalysis = "filtering")

      update_filter_ui(session)
      update_filter_combined_ui(session)
    }
    if (!is.null(rtdata$prediction$candidates)) { # Prediction

      rtdata$prediction$candidates <-
        create_other_isobars(rtdata$prediction$candidates,
                             typeanalysis = "prediction")

      update_prediction_ui(session)
    }

    # Notification that it has been imported successfully.
    showNotification(
      paste0("The object '", rtdata_name, "' was imported successfully."),
      type = "message",
      duration = NULL
    )

    # Remove the initial load notification.
    removeNotification(id_loading)

  })



  #============================================================================#
  # Load file with the classes and adducts to predict.
  #============================================================================#
  observeEvent(input$filedbs, {

    # Do nothing if default database is selected
    if (input$type_lipidclasses == "lipidclasses_default" || is.null(input$filedbs)) {
      return(NULL)
    }

    # Show loading message
    loading_id <- showNotification(
      HTML("Importing lipid class custom file... <span class='loader'></span>"),
      type = "default", duration = NULL
    )

    # Validate file input
    if (is.null(input$filedbs) || input$filedbs$datapath == "") {
      removeNotification(loading_id)
      showNotification("No file uploaded", type = "error")
      return(NULL)
    }


    tryCatch({

      filedbs_path <- input$filedbs$datapath
      filedbspathReactive(filedbs_path)

      # Initialize captured messages vector
      captured_messages <- character()

      generateadductsmzdbs <- withCallingHandlers({

        LipidRT::createadductsmzdbs(
          filedbs = filedbs_path,
          sep = input$sep_filedbs,
          adducts = adductsTableReactive()
        )

      }, message = function(m) {

        captured_messages <<- c(captured_messages, conditionMessage(m))
        invokeRestart("muffleMessage")

      })

      # Extract filedbs and dbs
      filedbs <- generateadductsmzdbs$filedbs
      dbs <- generateadductsmzdbs$adductsmzdbs

      View(filedbs)
      View(dbs)

      # Validate structure
      if (!is.list(dbs) || length(dbs) == 0) {
        removeNotification(loading_id)
        filedbsReactive(NULL)
        dbsReactive(NULL)

        showNotification(
          paste0(
            "Invalid format: the file must contain a 'classes' column, the 'adducts_pos' and/or 'adducts_neg' column according to the selected polarity, and the chain columns ",
            "'fas', 'sph', 'alkyl', and 'vinyl'. \n",
            "Please check the file structure and the selected separators."
          ),
          type = "error",
          duration = NULL
        )
      } else {
        filedbsReactive(filedbs)
        dbsReactive(dbs)

        if (!is.null(rtdataReactive())) {
          rtdata <- rtdataReactive()
          rtdata$metadata$prediction$filedbs <- filedbs
          rtdata$metadata$prediction$dbs <- dbs
          rtdataReactive(rtdata)
        }


        removeNotification(loading_id)

        # If messages were captured, show them in a styled banner
        if (length(captured_messages) > 0) {
          # Format each line according to its prefix
          # Split the message into individual lines
          lines <- unlist(strsplit(captured_messages, "\n"))

          # Format line by line
          full_message <- paste(
            sapply(lines, function(line) {

              # List lines starting with "-"
              if (grepl("^\\s*-", line)) {
                colon_pos <- regexpr(":", line)
                if (colon_pos > 0) {
                  left_part <- substr(line, 1, colon_pos)
                  left_part <- gsub("'", "", left_part)  # Remove single quotes
                  right_part <- substr(line, colon_pos + 1, nchar(line))
                  paste0(
                    "<span style='color: #c00000; font-weight:bold;'>", left_part, "</span>",
                    "<span style='color: black;'>", right_part, "</span>"
                  )
                } else {
                  line_clean <- gsub("'", "", line)
                  paste0("<span style='color: #c00000;'>", line_clean, "</span>")
                }

                # Lines starting with "Allowed"
              } else if (grepl("^Allowed", line)) {
                colon_pos <- regexpr(":", line)
                if (colon_pos > 0) {
                  left_part <- substr(line, 1, colon_pos)
                  right_part <- substr(line, colon_pos + 1, nchar(line))
                  paste0(
                    "<span style='color: #c00000; font-weight:bold;'>", left_part, "</span>",
                    "<span style='color: black;'>", right_part, "</span>"
                  )
                } else {
                  line
                }

                # Line starting with "Note:"
              } else if (grepl("^Note: Only valid adducts will be used. ", line)) {
                # Replace the final sentence
                # Replace everything after "Note:" with the new text
                line <- sub("Note: Only valid adducts will be used. .*$",
                            "Note: Only valid adducts will be used. To add additional valid adducts, load a custom adductsTable with the adducts of interest.",
                            line)
                paste0("<span style='color: #c00000; font-weight:bold;'>", line, "</span>")

                # Everything else black and bold
              } else {
                paste0("<span style='color: #c00000; font-weight:bold;'>", line, "</span>")
              }

            }),
            collapse = "<br>"
          )

          # Display as a warning banner
          showNotification(
            HTML(paste0(
              "<div style='white-space: pre-wrap;'>",
              full_message,
              "</div>"
            )),
            type = "warning",
            duration = NULL
          )

        } else {
          # Normal success notification
          showNotification(
            "Lipid class file imported successfully.",
            type = "message",
            duration = 5
          )
        }
      }

    }, error = function(e) {
      removeNotification(loading_id)
      filedbsReactive(NULL)
      dbsReactive(NULL)

      showNotification(
        paste("Error reading custom filedbs:", e$message),
        type = "error",
        duration = NULL
      )

      showNotification(
        "Tip: Check if the selected column delimiter and adduct or chain separators match the file format.",
        type = "warning",
        duration = NULL
      )
    })
  })



  observeEvent(
    list(input$sep_filedbs, adductsTableReactive()), {

      # Only if it is not the default database
      if (input$type_lipidclasses == "filedbs_default") return(NULL)

      # Only if a previous file path exists
      if (!is.null(filedbspathReactive())) {

        # Show loading message
        loading_id <- showNotification(
          HTML("Re-importing lipid class custom file... <span class='loader'></span>"),
          type = "default", duration = NULL
        )

        # Initialize vector to capture messages
        captured_messages <- character()

        tryCatch({

          filedbs_path <- filedbspathReactive()

          # Capture messages using withCallingHandlers
          generateadductsmzdbs <- withCallingHandlers({

            LipidRT::createadductsmzdbs(
              filedbs = filedbs_path,
              sep = input$sep_filedbs,
              adducts = adductsTableReactive()
            )

          }, message = function(m) {
            captured_messages <<- c(captured_messages, conditionMessage(m))
            invokeRestart("muffleMessage")
          })

          filedbs <- generateadductsmzdbs$filedbs
          dbs <- generateadductsmzdbs$adductsmzdbs

          # Validate structure
          if (!is.list(dbs) || length(dbs) == 0) {
            filedbsReactive(NULL)
            dbsReactive(NULL)
            removeNotification(loading_id)

            showNotification(
              paste0(
                "Invalid format: the file must contain a 'classes' column, the 'adducts_pos' and/or 'adducts_neg' column according to the selected polarity, and the chain columns ",
                "'fas', 'sph', 'alkyl', and 'vinyl'. \n",
                "Please check the file structure and the selected separators."
              ),
              type = "error",
              duration = NULL
            )

            # Show messages if any were captured
            if (length(captured_messages) > 0) {
              lines <- unlist(strsplit(captured_messages, "\n"))

              full_message <- paste(
                sapply(lines, function(line) {

                  # List lines starting with "-"
                  if (grepl("^\\s*-", line)) {
                    colon_pos <- regexpr(":", line)
                    if (colon_pos > 0) {
                      left_part <- substr(line, 1, colon_pos)
                      left_part <- gsub("'", "", left_part)
                      right_part <- substr(line, colon_pos + 1, nchar(line))
                      paste0(
                        "<span style='color: #c00000; font-weight:bold;'>", left_part, "</span>",
                        "<span style='color: black;'>", right_part, "</span>"
                      )
                    } else {
                      line_clean <- gsub("'", "", line)
                      paste0("<span style='color: #c00000;'>", line_clean, "</span>")
                    }

                    # Lines starting with "Allowed"
                  } else if (grepl("^Allowed", line)) {
                    colon_pos <- regexpr(":", line)
                    if (colon_pos > 0) {
                      left_part <- substr(line, 1, colon_pos)
                      right_part <- substr(line, colon_pos + 1, nchar(line))
                      paste0(
                        "<span style='color: #c00000; font-weight:bold;'>", left_part, "</span>",
                        "<span style='color: black;'>", right_part, "</span>"
                      )
                    } else {
                      line
                    }

                    # Line starting with "Note:"
                  } else if (grepl("^Note: Only valid adducts will be used. ", line)) {
                    line <- sub("Note: Only valid adducts will be used. .*$",
                                "Note: Only valid adducts will be used. To add additional valid adducts, load a custom adductsTable with the adducts of interest.",
                                line)
                    paste0("<span style='color: #c00000; font-weight:bold;'>", line, "</span>")

                    # Everything else black
                  } else {
                    paste0("<span style='color: black;'>", line, "</span>")
                  }

                }),
                collapse = "<br>"
              )

              showNotification(
                HTML(paste0("<div style='white-space: pre-wrap;'>", full_message, "</div>")),
                type = "warning",
                duration = NULL
              )
            }

          } else {
            # Successful import
            filedbsReactive(filedbs)
            dbsReactive(dbs)

            if (!is.null(rtdataReactive())) {
              rtdata <- rtdataReactive()
              rtdata$metadata$prediction$filedbs <- filedbs
              rtdata$metadata$prediction$dbs <- dbs
              rtdataReactive(rtdata)
            }

            removeNotification(loading_id)

            if (length(captured_messages) > 0) {
              # Show messages if any
              lines <- unlist(strsplit(captured_messages, "\n"))
              full_message <- paste(
                sapply(lines, function(line) {
                  if (grepl("^\\s*-", line)) {
                    colon_pos <- regexpr(":", line)
                    if (colon_pos > 0) {
                      left_part <- substr(line, 1, colon_pos)
                      left_part <- gsub("'", "", left_part)
                      right_part <- substr(line, colon_pos + 1, nchar(line))
                      paste0(
                        "<span style='color: #c00000; font-weight:bold;'>", left_part, "</span>",
                        "<span style='color: black;'>", right_part, "</span>"
                      )
                    } else {
                      line_clean <- gsub("'", "", line)
                      paste0("<span style='color: #c00000;'>", line_clean, "</span>")
                    }
                  } else if (grepl("^Allowed", line)) {
                    colon_pos <- regexpr(":", line)
                    if (colon_pos > 0) {
                      left_part <- substr(line, 1, colon_pos)
                      right_part <- substr(line, colon_pos + 1, nchar(line))
                      paste0(
                        "<span style='color: #c00000; font-weight:bold;'>", left_part, "</span>",
                        "<span style='color: black;'>", right_part, "</span>"
                      )
                    } else {
                      line
                    }
                  } else if (grepl("^Note: Only valid adducts will be used. ", line)) {
                    line <- sub("Note: Only valid adducts will be used. .*$",
                                "Note: Only valid adducts will be used. To add additional valid adducts, load a custom adductsTable with the adducts of interest.",
                                line)
                    paste0("<span style='color: #c00000; font-weight:bold;'>", line, "</span>")
                  } else {
                    paste0("<span style='color: black;'>", line, "</span>")
                  }
                }),
                collapse = "<br>"
              )
              showNotification(
                HTML(paste0("<div style='white-space: pre-wrap;'>", full_message, "</div>")),
                type = "warning",
                duration = NULL
              )
            } else {
              showNotification(
                "Lipid class file re-imported successfully using the new separators.",
                type = "message",
                duration = 5
              )
            }
          }

        }, error = function(e) {
          removeNotification(loading_id)
          filedbsReactive(NULL)
          dbsReactive(NULL)

          showNotification(
            paste("Error reading lipid class file:", e$message),
            type = "error",
            duration = NULL
          )

          showNotification(
            "Tip: Check if the selected column delimiter and adduct or chain separators match the file format.",
            type = "warning",
            duration = NULL
          )
        })
      }
    }
  )





  output$downloadfiledbs <- downloadHandler(
    filename = function(){paste(input$jobname, "filedbs.csv", sep="_")} ,
    content = function(file) {

      # Create the data frame using LipidRT
      filedbs_default <- LipidRT::defaultfiledbs(sep = input$sep_filedbs)

      # Save as CSV
      write.table(filedbs_default, file = file, sep = input$sep_filedbs,
                  row.names = FALSE, col.names = TRUE,  quote = TRUE )
    }
  )



  #============================================================================#
  # Import adductsTable from the global environment into the app.
  #============================================================================#
  observeEvent(input$fileadductsTable, {

    # Do nothing if default database is selected
    if (input$type_lipidclasses == "filedbs_default" || is.null(input$fileadductsTable)) {
      return(NULL)
    }

    # Show loading message
    loading_id <- showNotification(
      HTML("Importing adductsTable custom file... <span class='loader'></span>"),
      type = "default", duration = NULL
    )

    # Validate file input
    if (is.null(input$fileadductsTable) || input$fileadductsTable$datapath == "") {
      removeNotification(loading_id)
      showNotification("No file uploaded", type = "error")
      return(NULL)
    }


    tryCatch({

      fileadductsTablepath <- input$fileadductsTable$datapath
      adductsTablepathReactive(fileadductsTablepath)

      adductsTable_csv <- read.csv(fileadductsTablepath, header = TRUE, sep = input$sep_adductsTable,
                                   stringsAsFactors = FALSE, check.names = TRUE,
                                   na.strings = c("", "NA"))

      # If it has the correct structure...
      if (is.data.frame(adductsTable_csv) &&
          all(c("adduct", "mdiff", "charge", "n") %in% names(adductsTable_csv))) {

        adductsTableReactive(adductsTable_csv)

        if (!is.null(rtdataReactive())) {
        rtdata <- rtdataReactive()
        rtdata$metadata$prediction$adductsTable <- adductsTable_csv
        rtdataReactive(rtdata)
        }

        # Notification that it has been imported successfully.
        showNotification(
          paste0("The adductsTable file was imported successfully."),
          type = "message",
          duration = NULL
        )

        # Remove the initial load notification.
        removeNotification(loading_id)

      } else {
        # If it doesn’t have the correct structure, it notifies you and stops.
        showNotification(
          paste0("The adductsTable file does not have the required format."),
          type = "error"
        )
        removeNotification(loading_id)
        return()
      }


    }, error = function(e) {
      removeNotification(loading_id)

      showNotification(
        paste("Error reading custom adductsTable:", e$message),
        type = "error",
        duration = NULL
      )

      showNotification(
        "Tip: Check if the selected column delimiter match the file format.",
        type = "warning",
        duration = NULL
      )
    })
  })


  observeEvent(input$import_adductsTable_renv, {

    adductsTable_name <- input$name_adductsTable_renv  # Get the name of the loaded object.

    id_loading <- showNotification(
      HTML(paste0("Importing the object ", adductsTable_name, " from the global environment... <span class='loader'></span>")),
      type = "default",
      duration = NULL
    )

    # Check if the object exists in the global environment.
    if (!exists(adductsTable_name, envir = .GlobalEnv)) {
      showNotification(
        paste0("The object '", adductsTable_name, "' does not exist in the global environment."),
        type = "error",
        duration = NULL
      )
      removeNotification(id_loading)
      return()
    }

    # Get the object from the global environment.
    adductsTable_modif <- get(adductsTable_name, envir = .GlobalEnv)

    # If it has the correct structure...
    if (is.data.frame(adductsTable_modif) &&
        all(c("adduct", "mdiff", "charge", "n") %in% names(adductsTable_modif))) {

      adductsTableReactive(adductsTable_modif)


      if (!is.null(rtdataReactive())) {
      rtdata <- rtdataReactive()
      rtdata$metadata$prediction$adductsTable <- adductsTable_modif
      rtdataReactive(rtdata)
      }

    } else {
      # If it doesn’t have the correct structure, it notifies you and stops.
      showNotification(
        paste0("The object '", adductsTable_name, "' does not have the required format."),
        type = "error"
      )
      removeNotification(id_loading)
      return()
    }

    # Notification that it has been imported successfully.
    showNotification(
      paste0("The object '", adductsTable_name, "' was imported successfully."),
      type = "message",
      duration = NULL
    )

    # Remove the initial load notification.
    removeNotification(id_loading)

  })



  output$downloadadductsTable <- downloadHandler(
    filename = function(){paste(input$jobname, "adductsTable.csv", sep="_")} ,
    content = function(file) {

      # Create the data frame using LipidRT
      dbs <- LipidRT::assignDBs()
      adductsTable <- dbs$adductsTable

      # Save as CSV
      write.table(adductsTable, file = file, sep = input$sep_adductsTable,
                  row.names = FALSE, col.names = TRUE,  quote = TRUE )
    }
  )



  #============================================================================#
  # Export rtdata from the global environment to the app.
  #============================================================================#
  # Show modal when clicking the main button
  observeEvent(input$export_rtdata1, {
    showModal(modalDialog(
      title = "Export Object",
      textInput("export_object_name1", "Enter the name for the object:", value = "rtdata"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_export1", "Export")
      )
    ))
  })

  # Export button
  observeEvent(input$confirm_export1, {
    rtdata <- rtdataReactive()

    if (!is.null(rtdata$prediction$candidates$other_isobars)) {
      rtdata$prediction$candidates$other_isobars <- NULL
    }

    obj_name <- input$export_object_name1

    if (obj_name == "") {
      showNotification("Please enter a valid object name.", type = "error")
      return(NULL)
    }

    # Keep the global environment.
    envir <- ".GlobalEnv"
    if (envir == ".GlobalEnv") {
      envir <- tryCatch({
        get(envir, envir = .GlobalEnv)
      }, error = function(e) {
        NULL
      })
    } else {
      stop("Wrong envir. It must be '.GlobalEnv'.")
    }

    # Send the rtdata to the global environment.
    if (!is.null(rtdata)) {
        id_loading <- showNotification(HTML(paste0("Exporting the object '", obj_name, "' to the global environment... <span class='loader'></span>")),
                                       type = "default", duration = NULL)
        assign(obj_name, rtdata, envir = envir) # Default name "rtdata".

        removeNotification(id_loading)
        showNotification(paste0("The object '", obj_name, "' has been successfully exported."), type = "message", duration = NULL)
    } else {
      # If the object doesn’t exist, nothing will be exported.
        showNotification(paste0("The object '", obj_name, "'does not exist.", type = "error", duration = NULL))
    }
    removeModal() # Close the modal after exporting
  })


  observeEvent(input$export_rtdata2, {
    showModal(modalDialog(
      title = "Export Object",
      textInput("export_object_name2", "Enter the name for the object:", value = "rtdata"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_export2", "Export")
      )
    ))
  })

  # Export button
  observeEvent(input$confirm_export2, {
    rtdata <- rtdataReactive()

    if (!is.null(rtdata$prediction$candidates$other_isobars)) {
      rtdata$prediction$candidates$other_isobars <- NULL
    }

    obj_name <- input$export_object_name2

    if (obj_name == "") {
      showNotification("Please enter a valid object name.", type = "error")
      return(NULL)
    }

    # Keep the global environment.
    envir <- ".GlobalEnv"
    if (envir == ".GlobalEnv") {
      envir <- tryCatch({
        get(envir, envir = .GlobalEnv)
      }, error = function(e) {
        NULL
      })
    } else {
      stop("Wrong envir. It must be '.GlobalEnv'.")
    }

    # Send the rtdata to the global environment.
    if (!is.null(rtdata)) {
      id_loading <- showNotification(HTML(paste0("Exporting the object '", obj_name, "' to the global environment... <span class='loader'></span>")),
                                     type = "default", duration = NULL)
      assign(obj_name, rtdata, envir = envir) # Default name "rtdata".

      removeNotification(id_loading)
      showNotification(paste0("The object '", obj_name, "' has been successfully exported."), type = "message", duration = NULL)
    } else {
      # If the object doesn’t exist, nothing will be exported.
      showNotification(paste0("The object '", obj_name, "'does not exist.", type = "error", duration = NULL))
    }
    removeModal() # Close the modal after exporting
  })





  output$download_PeaklistResults_ui <- renderUI({
    rtdata <- rtdataReactive()

    # Only show if there is valid data
    has_valid_data <- !is.null(rtdata$features) && nrow(rtdata$features) > 0

    if (has_valid_data) {
      tags$div(
        downloadButton("downloadPeaklistResults", "Peaklist",
                       class = "btn-primary",
                       style = "color: #fff; background-color: #435f7a ;
                       border-color: #2c3e50 ; width: 100%;
                       text-align: center; white-space: normal;")
        )
    } else { NULL }
  })

  output$download_FilterResults_ui <- renderUI({
    rtdata <- rtdataReactive()

    # Only show if there is valid data
    has_valid_data <- !is.null(rtdata$filtering$data) &&
      length(rtdata$filtering$plots) > 0

    if (has_valid_data) {
      tags$div(
        downloadButton("downloadFilterResults", "Filter",
                       class = "btn-primary",
                       style = "color: #fff; background-color: #435f7a ;
                       border-color: #2c3e50 ; width: 100%;
                       text-align: center; white-space: normal;")
      )
    } else { NULL }
  })

  output$download_PredictionResults_ui <- renderUI({
    rtdata <- rtdataReactive()

    # Only show if there is valid data
    has_valid_data <- !is.null(rtdata$prediction$candidates) &&
      length(rtdata$prediction$plots) > 0

    if (has_valid_data) {
      tags$div(
        downloadButton("downloadPredictionResults", "Prediction",
                       class = "btn-primary",
                       style = "color: #fff; background-color: #435f7a ;
                       border-color: #2c3e50 ; width: 100%;
                       text-align: center; white-space: normal;")
      )
    } else { NULL }
  })


  output$download_TargetResults_ui <- renderUI({
    rtdata <- rtdataReactive()

    # Only show if there is valid data
    has_valid_data <- !is.null(rtdata$prediction$predictions) &&
      nrow(rtdata$prediction$predictions) > 0

    if (has_valid_data) {
      tags$div(
        downloadButton("downloadTargetResults", "Target Prediction",
                       class = "btn-primary",
                       style = "color: #fff; background-color: #435f7a ;
                       border-color: #2c3e50 ; width: 100%;
                       text-align: center; white-space: normal;")
      )
    } else { NULL }
  })

  # Make the download button appear only if data is loaded.
  output$download_rtdata1_ui <- renderUI({
    if (!is.null(rtdataReactive())) {
      downloadButton("download_rtdata1", ".RData", class = "btn-primary",
                     style = "color: #fff; background-color: #2c3e50 ;
                     border-color: #2c3e50 ; width: 100%;
                     text-align: center; white-space: normal;")
    }
  })

  # Make the export button appear only if data is loaded.
  output$export_rtdata1_ui <- renderUI({
    if (!is.null(rtdataReactive())) {
      actionButton("export_rtdata1",  "Export in R environment",
                   style = "color: #fff; background-color: #2c3e50 ;
                   border-color: #2c3e50 ; width: 100%;
                   text-align: center; white-space: normal;")
    }
  })

  # Make the download button appear only if data is loaded.
  output$download_rtdata2_ui <- renderUI({
    if (!is.null(rtdataReactive())) {
      downloadButton("download_rtdata2", ".RData", class = "btn-primary",
                     style = "color: #fff; background-color: #2c3e50 ;
                     border-color: #2c3e50 ; width: 100%;
                     text-align: center; white-space: normal;")
    }
  })

  # Make the export button appear only if data is loaded.
  output$export_rtdata2_ui <- renderUI({
    if (!is.null(rtdataReactive())) {
      actionButton("export_rtdata2",  "Export in R environment",
                   style = "color: #fff; background-color: #2c3e50 ;
                   border-color: #2c3e50 ; width: 100%;
                   text-align: center; white-space: normal;")
    }
  })


  #============================================================================#
  # Run Button
  #============================================================================#
  # observeEvent(input$Run_analysis, {
  #
  #   # Filtering/prediction check
  #   if ((input$apply_filter || input$apply_prediction) && is.null(rtdataReactive())) {
  #     showNotification(
  #       HTML("Cannot perform filtering or prediction: no data has been imported yet."),
  #       type = "error", duration = 6
  #     )
  #     return(NULL)
  #   }
  #
  #   if (!is.null(rtdataReactive())) {
  #     rtdata <- rtdataReactive()
  #
  #     ##########################################################################
  #     # Apply filtering ........................................................
  #     ##########################################################################
  #     if (input$apply_filter) {
  #       id_loading <- showNotification(HTML("Performing the filtering... <span class='loader'></span>"),
  #                                      type = "default", duration = NULL)
  #
  #       # Hide the filter_results section while filtering is running.
  #       output$filter_ui <- renderUI({ NULL })
  #       output$filter_combined_ui <- renderUI({ NULL })
  #
  #       lipid_classes_combined(NULL)
  #       selected_lipid_classes_combined(NULL)
  #
  #
  #       # Initialize captured messages vector
  #       captured_messages_filter <- character()
  #
  #       filter_output <- capture.output({
  #
  #         rtdata <- LipidRT::filterannotations(
  #           rtdata = rtdata,
  #           isomeric_plasmalogens = as.logical(input$isomeric_plasmalogens),
  #           filedbs = filedbspathReactive(),
  #           sep = input$sep_filedbs,
  #           adducts = adductsTableReactive(),
  #           dmz = as.numeric(input$dmz_predict),
  #           ppm = as.logical(input$ppm_predict),
  #           distC = as.numeric(input$distC_filter),
  #           minC_filter = as.numeric(input$minCfilter_filter),
  #           distribution = as.character(input$distribution_filter),
  #           rtscorethr = as.numeric(input$score_filter),
  #           rttol_id = as.numeric(input$rttol_id),
  #           rttol_name = as.numeric(input$rttol_name),
  #           filter_plot = TRUE,
  #           verbose = TRUE
  #         )
  #       })
  #
  #       rtdataReactive(rtdata)
  #
  #       # Update 'Filter Results' in the interface.
  #       update_features_ui(session)
  #       update_filter_ui(session) # Update the "Filter Results" graphs.
  #       update_filter_combined_ui(session) # Update the "Filter Results" graphs.
  #
  #       # Notification that the filtering has been completed.
  #       removeNotification(id_loading)
  #       showNotification("The filtering has been successfully completed.", type = "message", duration = NULL)
  #
  #       # If messages were captured, show them in a styled banner
  #       if (length(filter_output) > 0) {
  #
  #         # Remove ANSI codes
  #         filter_output_clean <- gsub("\033\\[[0-9;]*m", "", filter_output)
  #
  #         formatted_lines <- sapply(seq_along(filter_output_clean), function(i) {
  #
  #           line <- filter_output_clean[i]
  #
  #           # First line: larger font
  #           if (i == 1) {
  #             return(paste0(
  #               "<span style='font-size:16px;color:#669999;'><strong>",
  #               line,
  #               "</strong></span>"
  #             ))
  #           }
  #
  #           # Remaining lines: color the part before the first |
  #           if (grepl("\\|", line)) {
  #
  #             split_pos <- regexpr("\\|", line)
  #
  #             left_part  <- substr(line, 1, split_pos - 1)
  #             right_part <- substr(line, split_pos, nchar(line))
  #
  #             return(paste0(
  #               "<span style='color:black;'><strong>", left_part, "</strong></span>",
  #               right_part
  #             ))
  #           }
  #
  #           line
  #         })
  #
  #         full_message_filter <- paste(formatted_lines, collapse = "<br>")
  #
  #         showNotification(
  #           HTML(
  #             paste0(
  #               "<div style='white-space: pre; font-family: monospace;'>",
  #               full_message_filter,
  #               "</div>"
  #             )
  #           ),
  #           type = "warning",
  #           duration = NULL
  #         )
  #       }
  #
  #     }
  #
  #     ##########################################################################
  #     # Apply RT prediction ....................................................
  #     ##########################################################################
  #     if (input$apply_prediction) {
  #       id_loading <- showNotification(
  #         HTML("Performing the retention time prediction... <span class='loader'></span>"),
  #         type = "default", duration = NULL)
  #
  #       # Hide the prediction_results section while the prediction is running.
  #       output$prediction_ui <- renderUI({ NULL })
  #
  #       filedbs_path <- filedbspathReactive()
  #
  #       # Initialize captured messages vector
  #       captured_messages_prediction <- character()
  #
  #
  #       # print("Distribution: ")
  #       # print(as.character(input$distribution_predict))
  #       # print("subclasses: ")
  #       # print("subclasses: ")
  #       # print(input$lipidclasses_predict)
  #       # print("filedbs: ")
  #       # print(filedbs_path)
  #       # print("sep: ")
  #       # print(input$sep_filedbs)
  #       # print("adducts: ")
  #       # print(adductsTableReactive())
  #       # print("dmz: ")
  #       # print(as.numeric(input$dmz_predict))
  #       # print("ppm: ",as.logical(input$ppm_predict))
  #       # print(paste0("ppm: ",as.logical(input$ppm_predict)))
  #       # print("rtscorethr: ")
  #       # print(as.numeric(input$score_predict))
  #       # print("minC_predict: ")
  #       # print(as.numeric(input$minCpredict_predict))
  #       # print("rt_thr: ")
  #       # print(as.numeric(input$rt_thr_predict))
  #       # print("rttol_id: ")
  #       # print(as.numeric(input$rttol_id))
  #       # print("rttol_name: ")
  #       # print(as.numeric(input$rttol_name))
  #
  #
  #       prediction_output <- capture.output({
  #
  #         rtdata <-
  #           LipidRT::rtpredictannotations(rtdata = rtdata,
  #                                         distribution = as.character(input$distribution_predict),
  #                                         subclasses = input$lipidclasses_predict,
  #                                         filedbs = filedbs_path,
  #                                         sep = input$sep_filedbs,
  #                                         adducts = adductsTableReactive(),
  #                                         dmz = as.numeric(input$dmz_predict),
  #                                         ppm = as.logical(input$ppm_predict),
  #                                         rtscorethr = as.numeric(input$score_predict),
  #                                         minC_predict = as.numeric(input$minCpredict_predict),
  #                                         rt_thr = as.numeric(input$rt_thr_predict),
  #                                         rttol_id = as.numeric(input$rttol_id),
  #                                         rttol_name = as.numeric(input$rttol_name),
  #                                         predict_plot = TRUE,
  #                                         verbose = TRUE)
  #
  #       })
  #
  #       # rtdata$prediction$candidates <-
  #       #   create_other_isobars(rtdata$prediction$candidates)
  #
  #       rtdataReactive(rtdata)
  #
  #       # Update "Prediction Results" in the interface
  #       update_features_ui(session)
  #       update_prediction_ui(session) # Update the "Prediction Results" graphs
  #
  #       # Notification that the prediction has been completed.
  #       removeNotification(id_loading)
  #       showNotification("The retention time prediction has been successfully completed.", type = "message", duration = NULL)
  #
  #       # If messages were captured, show them in a styled banner
  #       if (length(prediction_output) > 0) {
  #
  #         # Remove ANSI codes
  #         prediction_output_clean <- gsub("\033\\[[0-9;]*m", "", prediction_output)
  #
  #         formatted_lines <- sapply(seq_along(prediction_output_clean), function(i) {
  #
  #           line <- prediction_output_clean[i]
  #
  #           # First line: larger font
  #           if (i == 1) {
  #             return(paste0(
  #               "<span style='font-size:16px;color:#669999;'><strong>",
  #               line,
  #               "</strong></span>"
  #             ))
  #           }
  #
  #           # Remaining lines: color the part before the first |
  #           if (grepl("\\|", line)) {
  #
  #             split_pos <- regexpr("\\|", line)
  #
  #             left_part  <- substr(line, 1, split_pos - 1)
  #             right_part <- substr(line, split_pos, nchar(line))
  #
  #             return(paste0(
  #               "<span style='color: black;'><strong>", left_part, "</strong></span>",
  #               right_part
  #             ))
  #           }
  #
  #           line
  #         })
  #
  #         full_message_prediction <- paste(formatted_lines, collapse = "<br>")
  #
  #         showNotification(HTML(paste0(
  #               "<div style='white-space: pre; font-family: monospace;'>",
  #               full_message_prediction, "</div>")),
  #           type = "warning", duration = NULL)
  #       }
  #
  #     }
  #   }
  #
  # })

  observeEvent(input$Run_filtering_analysis, {

    # Filtering/prediction check
    if (is.null(rtdataReactive())) {
      showNotification(
        HTML("Cannot perform filtering or prediction: no data has been imported yet."),
        type = "error", duration = 6
      )
      return(NULL)
    }

    if (!is.null(rtdataReactive())) {
      rtdata <- rtdataReactive()

      ##########################################################################
      # Apply filtering ........................................................
      ##########################################################################
      # if (input$apply_filter) {
        id_loading <- showNotification(HTML("Performing the filtering... <span class='loader'></span>"),
                                       type = "default", duration = NULL)

        # Hide the filter_results section while filtering is running.
        output$filter_ui <- renderUI({ NULL })
        output$filter_combined_ui <- renderUI({ NULL })

        lipid_classes_combined(NULL)
        selected_lipid_classes_combined(NULL)


        # Initialize captured messages vector
        captured_messages_filter <- character()

        filter_output <- capture.output({

          rtdata <- LipidRT::filterannotations(
            rtdata = rtdata,
            isomeric_plasmalogens = as.logical(input$isomeric_plasmalogens),
            filedbs = filedbspathReactive(),
            sep = input$sep_filedbs,
            adducts = adductsTableReactive(),
            dmz = as.numeric(input$dmz_predict),
            ppm = as.logical(input$ppm_predict),
            distC_filter = as.numeric(input$distC_filter),
            minSpecies_filter = as.numeric(input$minSpecies_filter),
            distribution = as.character(input$distribution_filter),
            rtscorethr = as.numeric(input$score_filter),
            rttol_id = as.numeric(input$rttol_id),
            rttol_name = as.numeric(input$rttol_name),
            filter_plot = TRUE,
            verbose = TRUE
          )
        })

        rtdataReactive(rtdata)

        # Update 'Filter Results' in the interface.
        update_features_ui(session)
        update_filter_ui(session) # Update the "Filter Results" graphs.
        update_filter_combined_ui(session) # Update the "Filter Results" graphs.

        # Notification that the filtering has been completed.
        removeNotification(id_loading)
        showNotification("The filtering has been successfully completed.", type = "message", duration = NULL)

        # If messages were captured, show them in a styled banner
        if (length(filter_output) > 0) {

          # Remove ANSI codes
          filter_output_clean <- gsub("\033\\[[0-9;]*m", "", filter_output)

          formatted_lines <- sapply(seq_along(filter_output_clean), function(i) {

            line <- filter_output_clean[i]

            # First line: larger font
            if (i == 1) {
              return(paste0(
                "<span style='font-size:16px;color:#669999;'><strong>",
                line,
                "</strong></span>"
              ))
            }

            # Remaining lines: color the part before the first |
            if (grepl("\\|", line)) {

              split_pos <- regexpr("\\|", line)

              left_part  <- substr(line, 1, split_pos - 1)
              right_part <- substr(line, split_pos, nchar(line))

              return(paste0(
                "<span style='color:black;'><strong>", left_part, "</strong></span>",
                right_part
              ))
            }

            line
          })

          full_message_filter <- paste(formatted_lines, collapse = "<br>")

          showNotification(
            HTML(
              paste0(
                "<div style='white-space: pre; font-family: monospace;'>",
                full_message_filter,
                "</div>"
              )
            ),
            type = "warning",
            duration = NULL
          )
        }

      # }
        updateNavbarPage(
          session,
          inputId = "inTabset",
          selected = "tab5"   # Results
        )

        updateTabsetPanel(
          session,
          inputId = "tabs",
          selected = "filter_results"
        )
    }

  })

  observeEvent(input$Run_prediction_analysis, {

    # Filtering/prediction check
    if (is.null(rtdataReactive())) {
      showNotification(
        HTML("Cannot perform filtering or prediction: no data has been imported yet."),
        type = "error", duration = 6
      )
      return(NULL)
    }

    if (!is.null(rtdataReactive())) {
      rtdata <- rtdataReactive()

      ##########################################################################
      # Apply RT prediction ....................................................
      ##########################################################################
      # if (input$apply_prediction) {
        id_loading <- showNotification(
          HTML("Performing the retention time prediction... <span class='loader'></span>"),
          type = "default", duration = NULL)

        # Hide the prediction_results section while the prediction is running.
        output$prediction_ui <- renderUI({ NULL })

        filedbs_path <- filedbspathReactive()

        # Initialize captured messages vector
        captured_messages_prediction <- character()


        prediction_output <- capture.output({

          rtdata <-
            LipidRT::rtpredictannotations(rtdata = rtdata,
                                          distribution = as.character(input$distribution_predict),
                                          subclasses = input$lipidclasses_predict,
                                          filedbs = filedbs_path,
                                          sep = input$sep_filedbs,
                                          adducts = adductsTableReactive(),
                                          dmz = as.numeric(input$dmz_predict),
                                          ppm = as.logical(input$ppm_predict),
                                          rtscorethr = as.numeric(input$score_predict),
                                          distC_predict = as.numeric(input$distC_predict),
                                          minSpecies_predict = as.numeric(input$minSpecies_predict),
                                          rt_thr = as.numeric(input$rt_thr_predict),
                                          rttol_id = as.numeric(input$rttol_id),
                                          rttol_name = as.numeric(input$rttol_name),
                                          predict_plot = TRUE,
                                          verbose = TRUE)

        })

        # rtdata$prediction$candidates <-
        #   create_other_isobars(rtdata$prediction$candidates)

        rtdataReactive(rtdata)

        # Update "Prediction Results" in the interface
        update_features_ui(session)
        update_prediction_ui(session) # Update the "Prediction Results" graphs

        # Notification that the prediction has been completed.
        removeNotification(id_loading)
        showNotification("The retention time prediction has been successfully completed.", type = "message", duration = NULL)

        # If messages were captured, show them in a styled banner
        if (length(prediction_output) > 0) {

          # Remove ANSI codes
          prediction_output_clean <- gsub("\033\\[[0-9;]*m", "", prediction_output)

          formatted_lines <- sapply(seq_along(prediction_output_clean), function(i) {

            line <- prediction_output_clean[i]

            # First line: larger font
            if (i == 1) {
              return(paste0(
                "<span style='font-size:16px;color:#669999;'><strong>",
                line,
                "</strong></span>"
              ))
            }

            # Remaining lines: color the part before the first |
            if (grepl("\\|", line)) {

              split_pos <- regexpr("\\|", line)

              left_part  <- substr(line, 1, split_pos - 1)
              right_part <- substr(line, split_pos, nchar(line))

              return(paste0(
                "<span style='color: black;'><strong>", left_part, "</strong></span>",
                right_part
              ))
            }

            line
          })

          full_message_prediction <- paste(formatted_lines, collapse = "<br>")

          showNotification(HTML(paste0(
            "<div style='white-space: pre; font-family: monospace;'>",
            full_message_prediction, "</div>")),
            type = "warning", duration = NULL)
        }

        updateNavbarPage(
          session,
          inputId = "inTabset",
          selected = "tab5"   # Results
        )

        updateTabsetPanel(
          session,
          inputId = "tabs",
          selected = "prediction_results"
        )
      # }
    }

  })


  #============================================================================#
  # -- "Peaklist Results" interface.
  #============================================================================#
  update_features_ui <- function(session) {
    rtdata <- rtdataReactive()

      # Render dynamic content in the UI.
        output$features_ui <- renderUI({
          tagList(
            tags$div(
              style = "margin-top: 20px;",
              DT::DTOutput("features_table")) # Table
          )
        })

      # Render the interactive table with the data from rtdata$features
        output$features_table <- DT::renderDT({
          DT::datatable(
            rtdata$features,
            options = list(
              pageLength = 15,
              scrollX = TRUE,
              dom = 'frtip',
              autoWidth = FALSE
            ),
            escape = FALSE,
            class = "display"
          )
        })
  }


  #============================================================================#
  # -- "Targeted Predictions Results" interface.
  #============================================================================#
  observeEvent(rtdataReactive(), {

    if (is.null(rtdataReactive()$prediction$predictions)) {
      targetpredictionsReactive(NULL)
      return(NULL)
    }

    predictions <- LipidRT::targetpredictions(predictions = rtdataReactive()$prediction$predictions)
    targetpredictionsReactive(predictions)
  })


  observeEvent(targetpredictionsReactive(), {
    update_targeted_ui(session)
  })


  update_targeted_ui <- function(session) {

    target_predictions <- targetpredictionsReactive()

    output$targeted_ui <- renderUI({
      req(target_predictions)

      tagList(
        tags$div(
          style = "margin-top: 20px; text-align:center;",
          DT::DTOutput("target_table")
        )
      )
    })

    output$target_table <- DT::renderDT({
      req(target_predictions)

      DT::datatable(
        target_predictions,
        options = list(
          pageLength = 20,
          scrollX = FALSE,
          dom = 'frtip',
          autoWidth = TRUE,
          columnDefs = list(
            list(className = 'dt-center', targets = "_all")
          )
        ),
        class = "display",
        escape = FALSE
      )
    })
  }


  #============================================================================#
  # -- "Filter Results" graphic interface.
  #============================================================================#
  update_filter_ui <- function(session) {

    rtdata <- rtdataReactive()

    # Extract DB values for each Lipid Class
    db_choices <- DBextraction(rtdata)
    dbChoicesReactive(db_choices) # The DBs of each lipid class
    dbSelectionsReactive(db_choices) # The DBs selected in each lipid class by the user


    # UI: Render content of the filter results tab .............................
    output$filter_ui <- renderUI({
      tagList(
        # Lipid class selection and graph download
        div(
          style = "margin-top: 0px; border-bottom: 1px solid #ddd; padding-bottom: 0px;",
          div(
            class = "no-width-change-tabs",
            fluidRow(
              column(width = 6,
                tagList(
                  tags$div(
                    style = "display: flex; justify-content: flex-end; gap: 10px; margin-top: 0px; margin-bottom: 0px;",

    tags$div(
      style = "display: flex; align-items: center; justify-content: center; height: 100%;",

      # Tutorial Image
      actionButton(
        inputId = "manual_filter_info_btn",
        label = tags$img(
          src = "iconohorizontal_buttoninfo.png",
          style = "height: auto; width: auto; margin-right: 10px; height: 50px;"
        ),
        style = "background-color: transparent; border: none; padding: 0;
        margin-right: 15px; cursor: pointer;"
      ),

      # Floating menu initially hidden
      tags$div(
        id = "manual_filter_dropdown_menu",
        style = " display: none; position: fixed; top: 50%; left: 50%;
        transform: translate(-50%, -50%); height: 90vh; width: 80vw;
        overflow-y: auto; background-color: rgba(255, 255, 255, 1);
        padding: 20px; font-size: 12px; text-align: right; border: none;
        border-radius: 5px; box-shadow: 0 4px 20px rgba(0,0,0,0.3); z-index: 1000;",

        # Zoom controls
        tags$div(
          style = "position: sticky; top: 0; background-color: rgba(255,255,255,0.95);
          z-index: 10; padding-bottom: 5px; text-align: right;",
          actionButton(
            "zoom_out_dropdown_filter", "-",
            style = "width: 25px; height: 25px; padding: 0; font-size: 16px; vertical-align: middle;"
          ),
          span(
            id = "zoom_level_dropdown_filter", "100%",
            style = "font-weight: bold; margin: 0 8px; vertical-align: middle;"
          ),
          actionButton(
            "zoom_in_dropdown_filter", "+",
            style = "width: 25px; height: 25px; padding: 0; font-size: 16px; vertical-align: middle;"
          )
        ),

        # Filtering tutorial image
        tags$img(
          id = "manualFilter_dropdown",
          src = "info_ManualFilter.png",
          style = "width: 90%; height: auto; display: block; margin: 0 auto;"
        ),

        # JavaScript script for zoom
        tags$script(HTML("
      var zoomDropdown_filter = 100;
      var imgDropdown_filter = document.getElementById('manualFilter_dropdown');
      var zoomLevelSpanDropdown_filter = document.getElementById('zoom_level_dropdown_filter');

      document.getElementById('zoom_in_dropdown_filter').onclick = function() {
        if (zoomDropdown_filter < 200) {
          zoomDropdown_filter += 10;
          imgDropdown_filter.style.width = zoomDropdown_filter + '%';
          zoomLevelSpanDropdown_filter.textContent = zoomDropdown_filter + '%';
        }
      };
      document.getElementById('zoom_out_dropdown_filter').onclick = function() {
        if (zoomDropdown_filter > 30) {
          zoomDropdown_filter -= 10;
          imgDropdown_filter.style.width = zoomDropdown_filter + '%';
          zoomLevelSpanDropdown_filter.textContent = zoomDropdown_filter + '%';
        }
      };
    "))
      ),

      # Script to show/hide the dropdown when clicking the circle
      tags$script(HTML("
          var manualFilterBtn = document.getElementById('manual_filter_info_btn');
          var manualFilterMenu = document.getElementById('manual_filter_dropdown_menu');

          manualFilterBtn.onclick = function(event) {
            event.stopPropagation();
            if (manualFilterMenu.style.display === 'none') {
              manualFilterMenu.style.display = 'block';
            } else {
              manualFilterMenu.style.display = 'none';
            }
          };

          document.addEventListener('click', function(event) {
            if (!manualFilterMenu.contains(event.target) && !manualFilterBtn.contains(event.target)) {
              manualFilterMenu.style.display = 'none';
            }
          });
        "))
          ),

                  )
                )
              )

            ),
            # Tabs for results and manual filtering
            tabsetPanel(
              type = "tabs",

              # Tab: Manual Filter
              tabPanel(
                "Manual Filter",
                div(
                  style = "margin-top: 20px; margin-bottom: 40px;",
                  fluidRow(
                    column(2,

                      # Range Score
                      div(class = "info-box",
                          div("Check the retained lipids with a score below the threshold. They may need further review.",
                              style = "text-align: center;",
                              class = "info-text"),
                          uiOutput("toggle_button_ui")
                      ),

                      tags$details(
                        style = "margin-bottom: 15px; padding: 8px 12px; border: 1px solid #dcdcdc;
                             border-radius: 6px; background-color: #fafafa;",

                        ## HEADER — todo clicable
                        tags$summary(
                          style = "cursor: pointer; outline: none;",

                          # Title
                          tags$div(
                            style = "font-size: 14px; font-weight: bold; color: #333; margin-bottom: 4px;",
                            icon("chart-area"),
                            HTML("&nbsp;Graph style")
                          ),

                          # Explanation (always visible)
                          tags$div(
                            style = "font-size: 11px; color: #555; margin-top: 2px; text-align: center; ",
                            "Choose which data layers are displayed."
                          ),

                          # Click hint (below explanation)
                          tags$div(
                            style = "text-align: center; font-size: 10px; color: #777; margin-top: 5px;",
                            "click to expand / collapse"
                          )
                        ),

                        ## CONTENT — solo el select
                        div(
                          style = "margin-top: 10px;",

                          selectInput(
                            inputId = "graph_style_filter",
                            label   = NULL,   # quitamos el título interno
                            choices = list(
                              "Kept Data + Removed Data" = c(
                                "Kept Data + Removed Data"
                              ),
                              "Single layers" = c(
                                "Kept Data",
                                "Removed Data"
                              )
                            ),
                            selected = "Kept Data + Removed Data",
                            width = "100%"
                          )
                        )
                      ),

                      # RT and C filtering (stacked, collapsible)
                      tags$details(
                        style = "margin-bottom: 15px; border: 1px solid #dcdcdc; border-radius: 6px;
                             background-color: #fafafa; padding: 8px 12px;",

                        ## Header
                        tags$summary(
                          style = "display: block; cursor: pointer; outline: none;",

                          # Línea superior: icono + título
                          tags$div(
                            style = "display: flex; align-items: center; font-size: 14px;
                                 font-weight: bold; color: #333; margin-bottom: 2px;",
                            icon("filter"),
                            HTML("&nbsp;&nbsp;RT - C Ranges")
                          ),

                          # Texto explicativo (siempre visible)
                          tags$div(
                            style = "font-size: 11px; color: #555; line-height: 1.3; margin-bottom: 2px; text-align: center;",
                            HTML(
                              "Use these controls to modify the Retention Time and Carbon Number ranges displayed."
                            )
                          ),

                          # Click hint CENTRADO debajo
                          tags$div(
                            style = "text-align: center; font-size: 10px; color: #777; margin-top: 5px;",
                            "click to expand / collapse"
                          )
                        ),

                        ## Content (STACKED)
                        div(
                          style = "margin-top: 10px;",

                          # RT range
                          div(
                            style = "padding: 4px; border: 1px dashed #d0d0d0; border-radius: 4px;
                                 background-color: #ffffff; margin-bottom: 10px;",
                            uiOutput("rt_range_filter_ui")
                          ),

                          # C range
                          div(
                            style = "padding: 4px; border: 1px dashed #d0d0d0; border-radius: 4px;
                                 background-color: #ffffff;",
                            uiOutput("c_range_filter_ui")
                          )
                        )
                      ),

                      # DB selector
                      tags$details(
                        id = "db_details_filter",
                        style = "margin-bottom: 20px; border: 1px solid #dcdcdc;
                             border-radius: 6px; background-color: #fafafa; padding: 8px 12px;",

                        tags$summary(
                          style = "cursor: pointer; outline: none;",

                          tags$div(
                            style = "font-size: 14px; font-weight: bold; color: #333; margin-bottom: 4px;",
                            icon("database"),
                            HTML("&nbsp;DB values")
                          ),

                          tags$div(
                            style = "font-size: 11px; color: #555; margin-top: 2px; text-align: center;",
                            HTML("Select the DB values affected by <b>Keep</b> and <b>Remove</b> actions.")
                          ),

                          tags$div(
                            style = "text-align: center; font-size: 10px; color: #777; margin-top: 5px;",
                            "click to expand / collapse"
                          )
                        ),

                        # SOLO el contenido interno es reactivo
                        div(
                          style = "margin-top: 15px;",
                          uiOutput("db_filter_content")
                        )
                      ),

                      div(
                        style = "margin-top: 15px; padding: 10px; border: 1px solid #dcdcdc;
                             border-radius: 6px; background-color: #fafafa;",

                        # Title
                        tags$div(
                          style = "font-size: 14px; font-weight: bold; margin-bottom: 4px;",
                          icon("sliders-h"),
                          HTML("&nbsp;Modify filtering")
                        ),

                        # Keep button
                        actionButton(
                          "keepFilterButton",
                          label = HTML('<span style="font-size:16px; margin-right:6px;">&#x25B2;</span><b style="font-size:14px;">Keep</b>'),
                          style = "width: 100%; margin-top: 6px; margin-bottom: 6px; background-color: #ecf0f1;
                               color: black; border: 1px solid #d0d0d0; padding: 4px 8px;"),

                        # Remove button
                        actionButton(
                          "removeFilterButton",
                          label = HTML('<span style=\"font-size:16px; margin-right:6px;\">&#x25B3;</span><b style=\"font-size:14px;\">Remove</b>'),
                          style = "width: 100%; margin-bottom: 8px; background-color: #ecf0f1;
                               color: black; border: 1px solid #d0d0d0; padding: 4px 8px;"
                        ),

                        # Save button
                        # actionButton(
                        #   "saveFilterButton",
                        #   label = HTML(paste(icon("save"), '<b style="font-size:14px; margin-left:6px;">Save</b>')),
                        #   style = "width: 100%; margin-bottom: 4px; background-color: #ecf0f1;
                        #        color: black; border: 1px solid #d0d0d0; padding: 4px 8px;"
                        # ),
                        uiOutput("saveFilterButtonUI"),

                        # Save hint text
                        tags$div(
                          style = "font-size: 11px; color: #666; margin-top: 2px; text-align: center;",
                          icon("info-circle"),
                          HTML("&nbsp;Changes are only applied after saving")
                        )
                      )

                    ),
                    column(8,

                        # Interactive Prediction Graph
                        div(
                          class = "full-height",
                          plotly::plotlyOutput("manualFilterPlot", height = "700px")
                        ),

                        ## Info box under the graph
                        # div(
                        #   style = "margin-top: 12px; padding: 10px 12px; border: 1px solid #dcdcdc;
                        # border-left: 4px solid #6c757d; border-radius: 4px; background-color: #f8f9fa;
                        # font-size: 12px; color: #333;",
                        #   icon("info-circle"),
                        #   tags$strong(" Filtering impact: "),
                        #   "If you modify the filtering, the predictions will be recalculated accordingly."
                        # )
                        uiOutput("filterImpactMessage_filter")

                    ),
                    column(2,
                      # Legend image of the filter graph. If you click on the image, it doubles in size. To return it to its original size, tap anywhere on the screen.
                      div(
                        class = "full-width",
                        tags$img(id = "zoomable-image_filter",src = "legend_manualFilter.png",
                          width = "90%",
                          style = "cursor: pointer; transition: transform 0.3s ease-in-out;
                          position: relative; top: 0; left: 0;"
                        ),
                        tags$script(HTML("
                            $(document).ready(function() {
                              $('#zoomable-image_filter').click(function(event) {
                                event.stopPropagation();
                                $(this).toggleClass('zoomed');
                              });
                              $(document).click(function() {
                                $('#zoomable-image_filter').removeClass('zoomed');
                              });
                            });
                          ")),
                                              tags$style(HTML("
                            #zoomable-image_filter {
                              max-width: 100%;
                            }
                            #zoomable-image_filter.zoomed {
                              transform: scale(2);
                              transform-origin: top left;
                              position: absolute;
                              top: 0;
                              left: 0;
                            }
                        "))
                      ),

                      div(
                        style = "margin-top: 20px; padding: 10px; border: 1px solid #dcdcdc;
                             border-radius: 6px; background-color: #fafafa;",

                        # Title
                        tags$div(
                          style = "font-size: 14px; font-weight: bold; margin-bottom: 4px;",
                          icon("highlighter"),
                          HTML("&nbsp;Highlight data")
                        ),

                        # Subtitle / explanation
                        tags$div(
                          style = "font-size: 11px; color: #555; margin-bottom: 8px; text-align: center;",
                          "Mark or Unmark lipids in the plot to highlight them."
                        ),

                        # Mark button
                        actionButton(
                          "markFilterButton",
                          label = HTML(paste(icon("tag"), '<b style="font-size: 13px; margin-left: 6px;">Mark</b>')),
                          style = "width: 100%; margin-bottom: 6px; background-color: #ecf0f1;
                               color: black; border: 1px solid #d0d0d0; padding: 4px 8px;"),

                        # Unmark button
                        actionButton(
                          "unmarkFilterButton",
                          label = HTML(paste(icon("times"),'<b style="font-size: 13px; margin-left: 6px;">Unmark</b>')),
                          style = "width: 100%; background-color: #ecf0f1; color: black;
                               border: 1px solid #d0d0d0; padding: 4px 8px; margin-bottom: 4px;"),

                        # Save mark button
                        actionButton(
                          "saveMarkFilterButton",
                          label = HTML(paste(icon("save"),'<b style="font-size: 13px; margin-left: 6px;">Save</b>')),
                          style = "width: 100%; background-color: #ecf0f1; color: black;
                               border: 1px solid #d0d0d0; padding: 4px 8px; margin-bottom: 4px;"),

                        # Save hint text
                        tags$div(
                          style = "font-size: 11px; color: #666; margin-top: 2px; text-align: center;",
                          icon("info-circle"),
                          HTML("&nbsp;Changes are only applied after saving")
                        )
                      )

                    )
                  )
                )
              ),

              # Tab: Filter Results
              tabPanel(
                "Filter Results",
                div(
                  style = "margin-top: 20px;",
                  plotOutput("filterPlot", height = "600px") # Automatic Filter graph
                )

              )
            )
          )
        )
      )
    })

    output$saveFilterButtonUI <- renderUI({

      pending_changes <- !is.null(selectedindices_filter()) &&
        length(selectedindices_filter()) > 0

      if (pending_changes) {

        actionButton(
          "saveFilterButton",
          label = HTML(
            paste(
              icon("save"),
              '<b style="font-size:14px; margin-left:6px;">Save</b>'
            )
          ),
          style = "
        width: 100%;
        margin-bottom: 4px;
        background-color: #d6ebff;
        color: #003b6f;
        border: 1px solid #9bc7ef;
        padding: 4px 8px;
        font-weight: bold;
      "
        )

      } else {

        actionButton(
          "saveFilterButton",
          label = HTML(
            paste(
              icon("save"),
              '<b style="font-size:14px; margin-left:6px;">Save</b>'
            )
          ),
          style = "
        width: 100%;
        margin-bottom: 4px;
        background-color: #ecf0f1;
        color: black;
        border: 1px solid #d0d0d0;
        padding: 4px 8px;
      "
        )

      }
    })


    output$filterImpactMessage_filter <- renderUI({

      rtdata <- rtdataReactive()

      show_message <-
        "Filter_final" %in% colnames(rtdata$filtering$data)

      if (show_message) {
        return(NULL)
      }

      div(
        style = "
      margin-top: 12px;
      padding: 10px 12px;
      border: 1px solid #dcdcdc;
      border-left: 4px solid #6c757d;
      border-radius: 4px;
      background-color: #f8f9fa;
      font-size: 12px;
      color: #333;
    ",
        icon("info-circle"),
        tags$strong(" Filtering impact: "),
        "If you modify the filtering, the final prediction check will be reset and predictions will need to be checked again."
      )
    })

    # Filtered DB Interface Values .............................................
    # output$db_filter_ui <- renderUI({
    #   req(dbSelectionsReactive())
    #   req(db_choices)
    #   req(input$lipid_class_selector)
    #
    #   db_selections <- dbSelectionsReactive()
    #   current_class <- input$lipid_class_selector
    #
    #   all_db_values      <- db_choices$filtering[[current_class]]
    #   selected_db_values <- db_selections$filtering[[current_class]]
    #
    #   # Robust numeric parsing + remove NAs
    #   choices_vec <- as.numeric(all_db_values)
    #   choices_vec <- sort(unique(choices_vec[!is.na(choices_vec)]))
    #
    #   tagList(
    #     checkboxGroupInput(
    #       "db_values_filter",
    #       "Select DB values:",
    #       choices  = choices_vec,
    #       selected = selected_db_values,
    #       inline   = TRUE
    #     ),
    #
    #     # Subtitle under the input when there are no choices
    #     if (length(choices_vec) == 0) {
    #       tags$small(
    #         style = "display:block; color:#777; margin-top:4px;",
    #         "No lipids to display"
    #       )
    #     }
    #   )
    # })


    output$db_filter_content <- renderUI({

      req(dbChoicesReactive(), dbSelectionsReactive(), input$lipid_class_selector)

      current_class <- input$lipid_class_selector

      choices_vec <- as.numeric(dbChoicesReactive()$filtering[[current_class]])
      choices_vec <- sort(unique(choices_vec[!is.na(choices_vec)]))

      selected_db_values <- intersect(
        dbSelectionsReactive()$filtering[[current_class]],
        choices_vec
      )

      div(
        checkboxGroupInput(
          "db_values_filter",
          label    = NULL,
          choices  = choices_vec,
          selected = selected_db_values,
          inline   = FALSE
        ),

        if (length(choices_vec) == 0) {
          tags$div(
            style = "font-size: 11px; color: #777; margin-top: 6px;",
            "No DB values available"
          )
        }
      )
    })


    # Observe user changes and update db_selections dynamically ................
    observeEvent(input$db_values_filter, {
      req(input$lipid_class_selector)

      db_selections <- dbSelectionsReactive()
      current_class <- input$lipid_class_selector

      # New user selection
      db_selections$filtering[[current_class]] <- input$db_values_filter
      # Save the user's new selection
      dbSelectionsReactive(db_selections)
    })


    # Render the automatic filter graphs .......................................
    output$filterPlot <- renderPlot({
      req(input$lipid_class_selector)

      lipid_class <- input$lipid_class_selector

      # Check if lipid_class is not NULL or NA
      req(!is.null(lipid_class), !is.na(lipid_class))

        if (!is.null(rtdata)) {
          # Check that the list of plots is not NULL
          plot_list <- rtdata$filtering$plots
          if (!is.null(plot_list) && !is.null(names(plot_list)) && lipid_class %in% names(plot_list)) {
            return(plot_list[[lipid_class]])
          } else {
            validate(need(FALSE, "No data available for the selected lipid class."))
          }
        }
    })


    # Filtered C Range Interface Values ........................................
    output$c_range_filter_ui <- renderUI({
      req(minC_filter())
      req(maxC_filter())
      c_min <- minC_filter()
      c_max <- maxC_filter()

      sliderInput("c_range_filter", "CN Range:", min = c_min, max = c_max,
                  value = c(c_min, c_max), step = 1, width = "80%")
    })


  # Filtered RT Range Interface Values .......................................
  output$rt_range_filter_ui <- renderUI({
    req(minRT_filter(), maxRT_filter())

    rt_min <- minRT_filter()
    rt_max <- maxRT_filter()

    div(
      style = "margin-bottom: 0px;",

      tags$div(
        style = "font-size: 12px; font-weight: bold; margin-bottom: 2px; text-align: center;",
        "Retention Time Range"
      ),

      sliderInput(
        "rt_range_filter",
        label = NULL,
        min = rt_min,
        max = rt_max,
        value = c(rt_min, rt_max),
        step = 10,
        width = "100%"
      )
    )
  })


  # Filtered C Range Interface Values ......................................
  output$c_range_filter_ui <- renderUI({
    req(minC_filter())
    req(maxC_filter())
    c_min <- minC_filter()
    c_max <- maxC_filter()

    div(
      style = "margin-bottom: 0px;",

      tagList(
        tags$div(
          style = "font-size: 12px; font-weight: bold; margin-bottom: 2px; text-align: center;",
          "Carbon Number Range"
        ),
        sliderInput(
          "c_range_filter",
          label = NULL,
          min = c_min,
          max = c_max,
          value = c(c_min, c_max),
          step = 1,
          width = "100%"
        )))
  })
}


  #============================================================================#
  # -- 'Manual Filter' Graph Controls
  #============================================================================#
  # REMOVE FILTER ..............................................................
  observeEvent(input$removeFilterButton, {

    # Modified lipid indices in the interactive filtering graph
    selected <- plotly::event_data("plotly_selected", source = "manualFilterPlot_source")

    rtdata <- rtdataReactive()

    # Keep the indices of selected lipids that match the user-selected DBs
    if (!is.null(selected)) {
      new_indices <- as.numeric(unique(selected$customdata))

      if (show_filterscore()) {
        indices <- intersect(
          new_indices,
          rtdata$filtering$data$index[
            rtdata$filtering$data$DB %in% input$db_values_filter &
              rtdata$filtering$data$Score_rtpred < rtdata$metadata$filtering$rtscorethr
          ]
        )
      } else {
        indices <- intersect(
          new_indices,
          rtdata$filtering$data$index[
            rtdata$filtering$data$DB %in% input$db_values_filter
          ]
        )
      }

      if ("Filter_final" %in% colnames(rtdata$filtering$data)) {

        rtdata$filtering$data$Filter_final[
          rtdata$filtering$data$index %in% indices
        ] <- "Removed manually by the user"

      } else {

        idx <- rtdata$filtering$data$index %in% indices

        rtdata$filtering$data$checkfilter[idx] <- ifelse(

          ## PRIORITY 1: Removed with RT modelling (keep it as is)
          rtdata$filtering$data$checkmodel[idx] == "Removed with RT modelling",

          "Removed with RT modelling",

          ifelse(

            ## PRIORITY 2: manual removal cases
            rtdata$filtering$data$Filter[idx] %in% c("Kept") |

              (rtdata$filtering$data$checkmodel[idx] %in% c("Recovered with RT modelling") &
                 rtdata$filtering$data$Filter[idx] %in% c(
                   "Confidence level among isobars",
                   "Removed by RT rules",
                   "In-source fragments",
                   "Removed distant data",
                   "Outlier",
                   "RT coherence between subclasses")),

            "Removed manually by the user",

            ## DEFAULT
            rtdata$filtering$data$Filter[idx]
          )
        )
      }

    }


    selected_indices <- selectedindices_filter()
    selected_indices <- c(selected_indices, indices)
    selectedindices_filter(selected_indices)

    # Update rtdataReactive
    rtdataReactive(rtdata)

    # Render the updated interactive graph
    output$manualFilterPlot <- plotly::renderPlotly({
      req(input$lipid_class_selector)
      renderFilterPlot()
    })
  })

  # KEEP FILTER ................................................................
  observeEvent(input$keepFilterButton, {

    # Modified lipid indices in the interactive filtering graph
    selected <- plotly::event_data("plotly_selected", source = "manualFilterPlot_source")

    rtdata <- rtdataReactive()

    # Keep the indices of selected lipids that match the user-selected DBs
    if (!is.null(selected)) {
      new_indices <- as.numeric(unique(selected$customdata))
      if (show_filterscore()) {
        indices <- intersect(
          new_indices,
          rtdata$filtering$data$index[
            rtdata$filtering$data$DB %in% input$db_values_filter &
              rtdata$filtering$data$Score_rtpred < rtdata$metadata$filtering$rtscorethr
          ]
        )
      } else {
        indices <- intersect(
          new_indices,
          rtdata$filtering$data$index[
            rtdata$filtering$data$DB %in% input$db_values_filter
          ]
        )
      }



      if ("Filter_final" %in% colnames(rtdata$filtering$data)) {

        rtdata$filtering$data$Filter_final[
          rtdata$filtering$data$index %in% indices
        ] <- "Kept manually by the user"

      } else {

        idx <- rtdata$filtering$data$index %in% indices

        checkmodel <- rtdata$filtering$data$checkmodel[idx]
        Filter     <- rtdata$filtering$data$Filter[idx]

        rtdata$filtering$data$checkfilter[idx] <- ifelse(

          ## PRIORITY 1: keep "Recovered with RT modelling"
          checkmodel == "Recovered with RT modelling",
          "Recovered with RT modelling",

          ifelse(

            ## PRIORITY 2: manual keep cases
            Filter %in% c("Confidence level among isobars", "Removed by RT rules",
                          "In-source fragments", "Removed distant data",
                          "Outlier", "RT coherence between subclasses") |
              (checkmodel == "Removed with RT modelling" & Filter == "Kept"),

            "Kept manually by the user",

            ## DEFAULT
            Filter
          )
        )
      }

    }

    selected_indices <- selectedindices_filter()
    selected_indices <- c(selected_indices, indices)
    selectedindices_filter(selected_indices)

    # Update rtdataReactive
    rtdataReactive(rtdata)

    # Render the updated interactive graph
    output$manualFilterPlot <- plotly::renderPlotly({
      req(input$lipid_class_selector)
      renderFilterPlot()
    })
  })


  # # MARK FILTER ................................................................
  observeEvent(input$markFilterButton, {

    showModal(modalDialog(
      title = NULL,

      div(
        style = "
        padding: 15px;
        background-color: #f8f9fa;
        border-radius: 8px;
      ",

        h4(icon("tag"), " Mark selected lipids"),
        p(
          style = "color: #666;",
          "Enter a comment that will be attached to all selected lipids."
        ),

        textAreaInput(
          "mark_text",
          label = NULL,
          placeholder = "Write your comment here...",
          rows = 3
        )
      ),

      footer = tagList(
        modalButton("Cancel"),
        actionButton(
          "confirm_mark_filter",
          "Apply comment",
          class = "btn-primary"
        )
      ),

      easyClose = TRUE
    ))
  })

  observe({
    shinyjs::toggleState(
      "confirm_mark_filter",
      condition = nzchar(trimws(input$mark_text))
    )
  })


  observeEvent(input$confirm_mark_filter, {

    removeModal()

    selected <- plotly::event_data(
      "plotly_selected",
      source = "manualFilterPlot_source"
    )

    rtdata <- rtdataReactive()

    if (!is.null(selected)) {

      new_indices <- as.numeric(unique(selected$customdata))

      if (show_filterscore()) {
        indices <- intersect(
          new_indices,
          rtdata$filtering$data$index[
            rtdata$filtering$data$DB %in% input$db_values_filter &
              rtdata$filtering$data$Score_rtpred < rtdata$metadata$filtering$rtscorethr
          ]
        )
      } else {
        indices <- intersect(
          new_indices,
          rtdata$filtering$data$index[
            rtdata$filtering$data$DB %in% input$db_values_filter
          ]
        )
      }

      # Usar el texto escrito por el usuario
      rtdata$filtering$data$comments[
        rtdata$filtering$data$index %in% indices
      ] <- input$mark_text
    }

    rtdataReactive(rtdata)

    output$manualFilterPlot <- plotly::renderPlotly({
      req(input$lipid_class_selector)
      renderFilterPlot()
    })
  })


  # UNMARK FILTER ..............................................................
  observeEvent(input$unmarkFilterButton, {

    # Modified lipid indices in the interactive prediction graph
    selected <- plotly::event_data("plotly_selected", source = "manualFilterPlot_source")

    rtdata <- rtdataReactive()

    # Keep the indices of selected lipids that match the user-selected DBs
    if (!is.null(selected)) {
      new_indices <- as.numeric(unique(selected$customdata))
      if (show_filterscore()) {
        indices <- intersect(
          new_indices,
          rtdata$filtering$data$index[
            rtdata$filtering$data$DB %in% input$db_values_filter &
              rtdata$filtering$data$Score_rtpred < rtdata$metadata$filtering$rtscorethr
          ]
        )
      } else {
        indices <- intersect(
          new_indices,
          rtdata$filtering$data$index[
            rtdata$filtering$data$DB %in% input$db_values_filter
          ]
        )
      }

      rtdata$filtering$data$comments[
        rtdata$filtering$data$index %in% indices
      ] <- ""

    }

    # Update 'rtdataReactive'
    rtdataReactive(rtdata)

    # Render the updated interactive graph
    output$manualFilterPlot <- plotly::renderPlotly({
      req(input$lipid_class_selector)
      renderFilterPlot()
    })
  })


  # SAVE MARK FILTER ................................................................
  observeEvent(input$saveMarkFilterButton, {
    notif_id <- showNotification(HTML("The changes from the Highlight are being saved...  <span class='loader'></span>"),
                                 type = "default", duration = NULL)

    # Save the changes in the reactive rtdata
    rtdata <- rtdataReactive()

    # Changes in filtering .....................................................
    rtdata <- LipidRT::organizertdata(rtdata = rtdata,
                                      typeanalysis = "filtering")

    # Redo the graphs with the changes
    rtdata <- LipidRT::filteringgraph(rtdata = rtdata)

    if (!is.null(rtdata$prediction$candidates)) {
      # Modify the feature matrix with the changes
      rtdata <- LipidRT::organizertdata(rtdata = rtdata,
                                        typeanalysis = "prediction")

      # Redo the graphs with the changes
      rtdata <- LipidRT::predictgraph(rtdata = rtdata,
                                      eval_predict_plot = TRUE)
    }

    # Update 'rtdataReactive' ..................................................
    rtdataReactive(rtdata)

    # Update feature table
    update_features_ui(session)

    # Render automatic graphs
    output$filterPlot <- renderPlot({
      req(input$lipid_class_selector)
      lipid_class <- input$lipid_class_selector

      rtdata$filtering$plots[[lipid_class]]
    })

    # Render the updated interactive graph
    output$manualFilterPlot <- plotly::renderPlotly({
      req(input$lipid_class_selector)

      renderFilterPlot()
    })

    # Notification: Saved successfully (lasts 3 seconds) .......................
    removeNotification(notif_id)  # Remove the 'saving' notification
    showNotification("The changes made in the Highlight have been successfully saved.",
                     type = "message", duration = 3)
  })





  # SAVE FILTER ................................................................
  observeEvent(input$saveFilterButton, {
    notif_id <- showNotification(HTML("The changes from the filtering are being saved...  <span class='loader'></span>"),
                                 type = "default", duration = NULL)

    # Save the changes in the reactive rtdata
    rtdata <- rtdataReactive()

    # Changes in filtering .....................................................
    # Rebuild models after saving the changes in filtering.

    if (!"Filter_final" %in% colnames(rtdata$filtering$data)) {
      rtdata <- LipidRT::checkremovedids(rtdata = rtdata,
                                         distribution = input$distribution_filter,
                                         rtscorethr = input$score_filter)
    }

    rtdata <- LipidRT::organizertdata(rtdata = rtdata,
                                      typeanalysis = "filtering")

    # Redo the graphs with the changes
    rtdata <- LipidRT::filteringgraph(rtdata = rtdata)


    rtdata$filtering$data <- create_other_isobars(rtdata$filtering$data,
                                                  typeanalysis = "filtering")

    # Update 'rtdataReactive' ..................................................
    rtdataReactive(rtdata)

    # Render automatic graphs
    output$filterPlot <- renderPlot({
      req(input$lipid_class_selector)
      lipid_class <- input$lipid_class_selector

      rtdata$filtering$plots[[lipid_class]]
    })

    if (!is.null(rtdata$prediction$candidates)) {
      rtdata <- LipidRT::organizertdata(rtdata = rtdata,
                                        typeanalysis = "prediction")

      # Redo the graphs with the changes
      rtdata <- LipidRT::predictgraph(rtdata = rtdata,
                                      eval_predict_plot = TRUE)

      rtdata$prediction$candidates <- create_other_isobars(rtdata$prediction$candidates,
                                                           typeanalysis = "prediction")

      # Render automatic graphs
      output$predictionPlot <- renderPlot({
        req(input$lipid_class_selector)
        lipid_class <- input$lipid_class_selector

        rtdata$prediction$plots[[lipid_class]]
      })

      # Render the updated interactive graph
      output$manualPredictionPlot <- plotly::renderPlotly({
        req(input$lipid_class_selector)

        renderPredictionGraph()
      })

    }

    # Render the updated interactive graph
    output$manualFilterPlot <- plotly::renderPlotly({
      req(input$lipid_class_selector)

      renderFilterPlot()
    })

    # Update 'rtdataReactive'
    rtdataReactive(rtdata)

    # Redo the predictions if they exist ........................................ (En principio esto no va a pasar porque al hacer las predicciones se hace el Filter final directamente, pero lo dejo por si acaso)
    if (!is.null(rtdata$prediction$candidates) &
        (!"Filter_final" %in% colnames(rtdata$filtering$data))) { # if rtdata$prediction exists...
      removeNotification(notif_id)
      notif_id <- showNotification(HTML("The changes from the filtering are being saved...  <br>
                                        Predictions are being recalculated for the classes affected by the filtering...<span class='loader'></span>"),
                                   type = "default", duration = NULL)

      # propose new rt candidates
      rtdata <- proposertcandidates(rtdata = rtdata,
                                    filedbs = filedbspathReactive(),
                                    sep = input$sep_filedbs,
                                    classes = input$lipidclasses_predict,
                                    adducts = adductsTableReactive(),
                                    dmz = as.numeric(input$dmz_predict),
                                    ppm = as.logical(input$ppm_predict))


      rtdata <- LipidRT::rtpredictannotationsdata(rtdata = rtdata,
                                                  distribution = as.character(input$distribution_predict),
                                                  classes = input$lipidclasses_predict,
                                                  rtscorethr = as.numeric(input$score_predict),
                                                  minC_predict = as.numeric(input$minCpredict_predict),
                                                  rt_thr = as.numeric(input$rt_thr_predict),
                                                  rttol_id = as.numeric(input$rttol_id),
                                                  rttol_name = as.numeric(input$rttol_name),
                                                  predict_plot = TRUE)

      rtdata <- LipidRT::organizertdata(rtdata = rtdata,
                                         typeanalysis = "prediction")

      # Update 'rtdataReactive'
      rtdataReactive(rtdata)

      # Update feature table
      update_features_ui(session)

      # Refresh the list of available DBs for filtering and prediction
      db_choices <- DBextraction(rtdata)
      dbChoicesReactive(db_choices)
      dbSelectionsReactive(db_choices)

      # Render the updated interactive graph
      output$manualPredictionPlot <- plotly::renderPlotly({
        renderPredictionGraph()
      })

      # Render graphs based on the selected LipidClas
      output$predictionPlot <- renderPlot({
        req(input$lipid_class_selector)
        lipid_class <- input$lipid_class_selector

        rtdata$prediction$plots[[lipid_class]]
      })

    }

    output$manualFilterPlot_combined <- plotly::renderPlotly({

      rtdata <- rtdataReactive()

      if (!is.null(rtdata$filtering$data)) {
        rtdata$filtering$data <- create_other_isobars(rtdata$filtering$data,
                                                      typeanalysis = "filtering")
      }
      if (!is.null(rtdata$prediction$candidates)) {
        rtdata$prediction$candidates <- create_other_isobars(rtdata$prediction$candidates,
                                                             typeanalysis = "prediction")
      }

      rtdataReactive(rtdata)

      req(lipid_classes_combined())
      req(selected_lipid_classes_combined())

      lipid_classes_combined()
      selected_lipid_classes_combined()

      renderFilterPlot_combined()
    })

    selectedindices_filter(NULL)

    # Update feature table
    update_features_ui(session)

    # Notification: Saved successfully (lasts 3 seconds) .......................
    removeNotification(notif_id)  # Remove the 'saving' notification
    showNotification("The changes made in the filtering have been successfully saved.",
                     type = "message", duration = 3)

    # if (!"Filter_final" %in% colnames(rtdata$filtering$data)) {
    # showNotification("The predictions have been recalculated based on the updated filtering criteria.",
    #                  type = "message", duration = 3)
    # }
  })



  #============================================================================#
  # -- Filtered Interactive Graph
  #============================================================================#
  renderFilterPlot <- function() {
    req(input$lipid_class_selector)


    rtdata <- rtdataReactive()

    req(input$lipid_class_selector)
    lipidclass <- input$lipid_class_selector
    score <- rtdata$metadata$filtering$rtscorethr

    # Filtering data ...........................................................
    if ("Filter_final" %in% colnames(rtdata$filtering$data)) {
      rtdata$filtering$data$checkfilter <- rtdata$filtering$data$Filter_final
    }

    filter_data <- reactive({

      rtdata$filtering$data[
        rtdata$filtering$data$LipidSubclass == lipidclass &
          !(rtdata$filtering$data$checkfilter %in% c("Ignore duplicated",
                                                     "Not pass QC", "Unsupported class")),
      ]
    })

    # Filter data for the graph
    shape_levels <- c("Kept" = 16,
                      "Kept manually by the user" = 17,
                      "Recovered with RT modelling" = 15,
                      "Removed with RT modelling" = 0,
                      "In-source fragments" = 8,
                      "Removed distant data" = 5,
                      "Confidence level among isobars" = 4,
                      "Removed by RT rules" = 1,
                      "Removed manually by the user" = 2,
                      "Outlier" = 6,
                      "RT coherence between subclasses" = 3)

    filter_data_plot <- filter_data()[
      filter_data()$checkfilter %in% c("Kept",
                                       "Kept manually by the user",
                                       "Recovered with RT modelling",
                                       "Removed with RT modelling",
                                       "In-source fragments",
                                       "Removed distant data",
                                       "Confidence level among isobars",
                                       "Removed by RT rules",
                                       "Removed manually by the user",
                                       "Outlier",
                                       "RT coherence between subclasses"), , drop = FALSE]

    filter_data_plot$checkfilter <- factor(filter_data_plot$checkfilter, levels = names(shape_levels),
                                           exclude = NULL)

    req(nrow(filter_data_plot) > 0)
    req(!is.null(input$c_range_filter))
    req(!is.null(input$rt_range_filter))
    req(length(input$c_range_filter) == 2)
    req(length(input$rt_range_filter) == 2)


    # Set the default alpha value to 1 for all DBs .............................
    if (show_filterscore()) {
      filter_data_plot$alpha <- ifelse(filter_data_plot$DB %in% input$db_values_filter &
                                         filter_data_plot$Score_rtpred < score &
                                         filter_data_plot$Filter %in% c("Kept"), 1, 0.8)
    } else {
      filter_data_plot$alpha <- ifelse(filter_data_plot$DB %in% input$db_values_filter, 1, 0.8)
    }


    # Range of minimum and maximum values of C from the data. ..................
    # Range of C safely
    if (nrow(filter_data_plot) == 0) return(NULL)

    min_c <- min(filter_data_plot$C, na.rm = TRUE) - 1
    max_c <- max(filter_data_plot$C, na.rm = TRUE) + 1

    if (!is.finite(min_c) || !is.finite(max_c)) return(NULL)

    min_c <- ifelse(min_c %% 2 == 0, min_c - 2, min_c - 1)
    max_c <- ifelse(max_c %% 2 == 0, max_c + 2, max_c + 1)


    minC_filter(min_c)
    maxC_filter(max_c)

    minc_plot <- as.numeric(input$c_range_filter[1])
    maxc_plot <- as.numeric(input$c_range_filter[2])

    req(is.finite(minc_plot), is.finite(maxc_plot))
    req(length(minc_plot) == 1, length(maxc_plot) == 1)


    # Range of minimum and maximum RT values from the data .....................
    min_rt <- floor(min(filter_data_plot$RTsec, na.rm = TRUE))
    max_rt <- ceiling(max(filter_data_plot$RTsec, na.rm = TRUE))

    if (!is.finite(min_rt) || !is.finite(max_rt)) return(NULL)


    minRT_filter(min_rt)
    maxRT_filter(max_rt)

    minrt_plot <- as.numeric(input$rt_range_filter[1])
    maxrt_plot <- as.numeric(input$rt_range_filter[2])

    req(is.finite(minrt_plot), is.finite(maxrt_plot))
    req(length(minrt_plot) == 1, length(maxrt_plot) == 1)


    # selected graph style .....................................................
    if (input$graph_style_filter == "Kept Data + Removed Data") {
      filter_data_plot <-
        filter_data_plot[!(filter_data_plot$checkfilter %in% c("Not pass QC",
                                                             "Unsupported class")), ]
    } else if (input$graph_style_filter == "Kept Data") {
      filter_data_plot <-
        filter_data_plot[filter_data_plot$checkfilter %in% c("Kept",
                                                             "Kept manually by the user",
                                                             "Recovered with RT modelling"), ]
    } else if (input$graph_style_filter == "Removed Data") {
      filter_data_plot <-
        filter_data_plot[filter_data_plot$checkfilter %in% c("Removed with RT modelling",
                                                             "In-source fragments",
                                                             "Removed distant data",
                                                             "Confidence level among isobars",
                                                             "Removed by RT rules",
                                                             "Removed manually by the user",
                                                             "Outlier",
                                                             "RT coherence between subclasses"), ]
    }

    if (is.null(filter_data_plot$comments)) {
      filter_data_plot$comments <- ""
    }


    # Reactive to update the predictions .......................................
    # updated_predictions <- reactive({
    #
    #   req(rtdataReactive())
    #   rtdata <- rtdataReactive()
    #
    #   rtdata$filtering$predictions[
    #     rtdata$filtering$predictions$LipidSubclass == lipidclass,
    #   ]
    #
    # })

    # prediction lines .........................................................
    # check if prediction layer is required
    prediction_layer <- LipidRT::predictionlayergraph(rtdata = rtdata,
                                             lipidsubclass = lipidclass,
                                             min_c = minc_plot,
                                             max_c = maxc_plot,
                                             min_rt = minrt_plot,
                                             max_rt = maxrt_plot,
                                             alpha = 0.6,
                                             typeanalysis = "filtering")


    # otrer isobars
    # Split por LipidRTid
    # filter_data_plot$other_isobars_text <- vapply(
    #   filter_data_plot$other_isobars,
    #   function(x) {
    #
    #     if (is.na(x) || x == "") {
    #       return("")
    #     }
    #
    #     parts <- strsplit(x, "\\|")[[1]]
    #
    #     lines <- vapply(parts, function(y) {
    #       p <- strsplit(y, ",")[[1]]
    #       paste0(
    #         p[1], " --- ", p[2], "<br>",
    #         " &nbsp;&nbsp; - RTsec: ", p[3], "<br>",
    #         " &nbsp;&nbsp; - RTsec_rtpred: ", p[4], "<br>",
    #         " &nbsp;&nbsp; - Score_rtpred: ", p[5], "<br>"
    #       )
    #     }, character(1))
    #
    #     paste0(
    #       "<br><b>- Other isobars:</b>",
    #       "<br>   ",
    #       paste(lines, collapse = "<br>   ")
    #     )
    #   },
    #   character(1)
    # )
    filter_data_plot$other_isobars_text <- vapply(
      filter_data_plot$other_isobars,
      function(x) {

        if (is.na(x) || x == "") {
          return("")
        }

        parts <- strsplit(x, "\\|")[[1]]

        lines <- vapply(seq_along(parts), function(i) {
          p <- strsplit(parts[i], ",")[[1]]
          paste0(
            "<b> Other Isobar ", i, " -------------------- ", "</b><br>",
            " - ",p[1], " --- ", p[2], "<br>",
            " - RTsec_rtpred: ", p[4], "<br>",
            " - Score_rtpred: ", p[5], "<br>",
            "<b> ----------------------------------------- ", "</b><br>"
          )
        }, character(1))

        paste0(
          "<br>",
          paste(lines, collapse = "<br>")
        )
      },
      character(1)
    )


    # otrer isobars
    # Split por LipidRTid
    # filter_data_plot$other_isobars_text <- vapply(
    #   filter_data_plot$other_isobars,
    #   function(x) {
    #
    #     if (is.na(x) || x == "") {
    #       return("")
    #     }
    #
    #     parts <- strsplit(x, "\\|")[[1]]
    #
    #     lines <- vapply(parts, function(y) {
    #       p <- strsplit(y, ",")[[1]]
    #       paste0(
    #         p[1], " --- ", p[2],
    #         " (", p[3], ")"
    #       )
    #     }, character(1))
    #
    #     paste0(
    #       "<br><b>- Other isobars:</b>",
    #       "<br>   ",
    #       paste(lines, collapse = "<br>   ")
    #     )
    #   },
    #   character(1)
    # )


    # Original Label
    filter_data_plot$original_label_text <- mapply(
      function(filter_value, checkfilter_value) {

        if (!checkfilter_value %in% c(
          "Kept manually by the user",
          "Removed manually by the user",
          "Recovered with RT modelling",
          "Removed with RT modelling"
        )) {
          return("")
        }

        if (is.na(filter_value) || filter_value == "") {
          return("")
        }

        paste0(
          "<br><b>- Original label:</b> ",
          filter_value
        )
      },
      filter_data_plot$Filter,
      filter_data_plot$checkfilter,
      USE.NAMES = FALSE
    )

    filter_data_plot$adduct_text <- vapply(
      filter_data_plot$adducts_expanded,
      function(x) {

        if (is.na(x) || x == "") {
          return("")
        }

        adducts <- strsplit(x, "\\|")[[1]]

        if (length(adducts) == 1) {
          return(
            paste0("<b>", adducts[1], "</b>")
          )
        }

        adducts[1] <- paste0("<b>", adducts[1], "</b>")

        paste(adducts, collapse = " | ")
      },
      character(1)
    )

    # Interactive graphics .....................................................
    lipidname <- IS <- C <- mz <- RTsec <- DB <- index <- LipidRTid <- RTsec_rtpred <- confidenceLevel <- checkfilter <- c()
    # If all the DBs are selected...
    if(length(setdiff(sort(unique(filter_data()$DB)), sort(input$db_values_filter))) == 0 &
       !show_filterscore()) {

      p <- ggplot(data = filter_data_plot,
                  aes(x = C, y = RTsec, color = factor(DB), shape = factor(checkfilter),
                      text = paste0(
                        "<br><b>", lipidname, "  ---  </b>", adduct_text,
                        ifelse(IS == "IS","<br><b> - Internal Standard</b>",""),

                        # Line
                        "<br>",
                        paste0(
                          rep("─", max(1, nchar(
                            ifelse(IS == "IS",
                                   " - Internal Standard",
                                   paste0(ifelse(is.na(lipidname), "", lipidname), "  ", ifelse(is.na(adduct), "", adduct)))
                          ))),
                          collapse = ""
                        ),
                        "<br><b> - LipidRTid:</b> ", LipidRTid,
                        "<br><b> - mz:</b> ", round(mz,3),
                        "<br><b> - RTsec:</b> ", round(RTsec, 3),
                        "<br><b> - RTsec_rtpred:</b> ", round(RTsec_rtpred, 3),
                        "<br><b> - Score_rtpred:</b> ", Score_rtpred,
                        "<br><b> - confidenceLevel:</b> ", confidenceLevel,
                        "<br><b> - Filter:</b> ", checkfilter,

                        # original label
                        ifelse(
                          !is.na(original_label_text) & original_label_text != "",
                          paste0(
                            "<br>",
                            paste0(
                              rep("─", max(1, nchar(
                                ifelse(IS == "IS",
                                       " - Internal Standard",
                                       paste0(ifelse(is.na(lipidname), "", lipidname),
                                              "  ",
                                              ifelse(is.na(adduct), "", adduct)))
                              ))),
                              collapse = ""
                            ),
                            original_label_text
                          ),
                          ""
                        ),

                        # other isobars
                        ifelse(
                          !is.na(other_isobars_text) & other_isobars_text != "",
                          paste0(
                            "<br>",
                            paste0(
                              rep("─", max(1, nchar(
                                ifelse(IS == "IS",
                                       " - Internal Standard",
                                       paste0(ifelse(is.na(lipidname), "", lipidname),
                                              "  ",
                                              ifelse(is.na(adduct), "", adduct)))
                              ))),
                              collapse = ""
                            ),
                            other_isobars_text
                          ),
                          ""
                        ),

                        # highlight
                        ifelse(
                          !is.na(comments) & comments != "",
                          paste0(
                            "<br>",
                            paste0(
                              rep("─", max(1, nchar(
                                ifelse(IS == "IS",
                                       " - Internal Standard",
                                       paste0(ifelse(is.na(lipidname), "", lipidname),
                                              "  ",
                                              ifelse(is.na(adduct), "", adduct)))
                              ))),
                              collapse = ""
                            ),
                            "<br> [ <b>!</b> ] <b>", comments, "</b> [ <b>!</b> ] "
                          ),
                          ""
                        )
                        ),

                      customdata = index)) +
        geom_point(size = 2.2, stroke = 0.6) +
        labs(x = "Carbon Number", y = "Retention time (s)", shape = "DB") +
        theme_minimal() +
        theme(legend.position = "none") +
        scale_color_manual(name = "Double \nbonds",
                           values = c("#212121", "#EE2C2C", "#FF9500", "#FFE100",
                                      "#90EE90", "#228B22", "#B0E2FF", "#00B2EE",
                                      "#3A5FCD", "#BF3EFF", "#68228B", "#795548",
                                      "#9E9E9E", "#FFB6C1", "#FFD39B", "#FFF68F",
                                      "#B8860B", "#8B864E", "#B87F9E", "#7FFF00",
                                      "#FF69B4", "#DEBF4B", "#51FFFC", "#EAC4FF",
                                      "#D9FF68",
                                      "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
                                      "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
                                      "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173",
                                      "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
                                      "#969696", "#DD1C77", "#6BAED6", "#74C476", "#9E9AC8",
                                      "#E7969C", "#CEDB9C", "#9C9EDE", "#DBDB8D", "#C49C94",
                                      "#F7B6D2", "#C7C7C7", "#AEC7E8", "#FFBB78", "#98DF8A",
                                      "#FF9896", "#C5B0D5", "#C49C94", "#FDD0A2", "#BDBDBD",
                                      "#5254A3", "#6B6ECF", "#9C9EDE", "#637939", "#B5CF6B",
                                      "#CEDB9C", "#8CA252", "#BD9E39", "#E7BA52", "#E7CB94",
                                      "#843C39", "#AD494A", "#D6616B", "#E7969C", "#7B4173",
                                      "#A55194", "#CE6DBD", "#DE9ED6", "#3182BD", "#6BAED6",
                                      "#9ECAE1", "#C6DBEF", "#E6550D", "#FD8D3C", "#FDAE6B",
                                      "#FDD0A2", "#31A354", "#74C476", "#A1D99B", "#C7E9C0",
                                      "#756BB1", "#9E9AC8", "#BCBDDC", "#DADAEB", "#636363",
                                      "#969696", "#BDBDBD", "#D9D9D9", "#393B79", "#5254A3",
                                      "#6B6ECF", "#9C9EDE", "#637939", "#8CA252", "#B5CF6B",
                                      "#CEDB9C", "#8C6D31", "#BD9E39", "#E7BA52", "#E7CB94"),
                           breaks = seq(0:24)-1) +
        scale_shape_manual(values = shape_levels) +
        scale_x_continuous(limits = c(minc_plot, maxc_plot), breaks = seq(minc_plot, maxc_plot, by = 2)) +
        scale_y_continuous(limits = c(minrt_plot, maxrt_plot)) +

      geom_point(
        data = subset(filter_data_plot, !is.na(comments) & comments != ""),
        aes(x = C, y = RTsec),
        shape = 22,
        size = 5,
        fill = "black",
        alpha = 0.25,
        color = NA,
        inherit.aes = FALSE
      )

    } else { # If any of the DBs are not selected, take alpha into account
      p <- ggplot(data = filter_data_plot,
                  aes(x = C, y = RTsec, color = factor(DB), shape = factor(checkfilter),
                      text = paste0(
                        "<br><b>", lipidname, "  ---  </b>", adduct_text,
                        ifelse(IS == "IS","<br><b> - Internal Standard</b>",""),

                        # Line
                        "<br>",
                        paste0(
                          rep("─", max(1, nchar(
                            ifelse(IS == "IS",
                                   " - Internal Standard",
                                   paste0(ifelse(is.na(lipidname), "", lipidname), "  ", ifelse(is.na(adduct), "", adduct)))
                          ))),
                          collapse = ""
                        ),
                        "<br><b> - LipidRTid:</b> ", LipidRTid,
                        "<br><b> - mz:</b> ", round(mz,3),
                        "<br><b> - RTsec:</b> ", round(RTsec, 3),
                        "<br><b> - RTsec_rtpred:</b> ", round(RTsec_rtpred, 3),
                        "<br><b> - Score_rtpred:</b> ", Score_rtpred,
                        "<br><b> - confidenceLevel:</b> ", confidenceLevel,
                        "<br><b> - Filter:</b> ", checkfilter,

                        # original label
                        ifelse(
                          !is.na(original_label_text) & original_label_text != "",
                          paste0(
                            "<br>",
                            paste0(
                              rep("─", max(1, nchar(
                                ifelse(IS == "IS",
                                       " - Internal Standard",
                                       paste0(ifelse(is.na(lipidname), "", lipidname),
                                              "  ",
                                              ifelse(is.na(adduct), "", adduct)))
                              ))),
                              collapse = ""
                            ),
                            original_label_text
                          ),
                          ""
                        ),

                        # other isobars
                        ifelse(
                          !is.na(other_isobars_text) & other_isobars_text != "",
                          paste0(
                            "<br>",
                            paste0(
                              rep("─", max(1, nchar(
                                ifelse(IS == "IS",
                                       " - Internal Standard",
                                       paste0(ifelse(is.na(lipidname), "", lipidname),
                                              "  ",
                                              ifelse(is.na(adduct), "", adduct)))
                              ))),
                              collapse = ""
                            ),
                            other_isobars_text
                          ),
                          ""
                        ),

                        # highlight
                        ifelse(
                          !is.na(comments) & comments != "",
                          paste0(
                            "<br>",
                            paste0(
                              rep("─", max(1, nchar(
                                ifelse(IS == "IS",
                                       " - Internal Standard",
                                       paste0(ifelse(is.na(lipidname), "", lipidname),
                                              "  ",
                                              ifelse(is.na(adduct), "", adduct)))
                              ))),
                              collapse = ""
                            ),
                            "<br> [ <b>!</b> ] <b>", comments, "</b> [ <b>!</b> ] "
                          ),
                          ""
                        )
                        ),

                      alpha = alpha, customdata = index)) +
        geom_point(size = 2.2, stroke = 0.6) +
        labs(x = "Carbon Number", y = "Retention time (s)", shape = "DB") +
        theme_minimal() +
        theme(legend.position = "none") +
        scale_color_manual(name = "Double \nbonds",
                           values = c("#212121", "#EE2C2C", "#FF9500", "#FFE100",
                                      "#90EE90", "#228B22", "#B0E2FF", "#00B2EE",
                                      "#3A5FCD", "#BF3EFF", "#68228B", "#795548",
                                      "#9E9E9E", "#FFB6C1", "#FFD39B", "#FFF68F",
                                      "#B8860B", "#8B864E", "#B87F9E", "#7FFF00",
                                      "#FF69B4", "#DEBF4B", "#51FFFC", "#EAC4FF",
                                      "#D9FF68",
                                      "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
                                      "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
                                      "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173",
                                      "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
                                      "#969696", "#DD1C77", "#6BAED6", "#74C476", "#9E9AC8",
                                      "#E7969C", "#CEDB9C", "#9C9EDE", "#DBDB8D", "#C49C94",
                                      "#F7B6D2", "#C7C7C7", "#AEC7E8", "#FFBB78", "#98DF8A",
                                      "#FF9896", "#C5B0D5", "#C49C94", "#FDD0A2", "#BDBDBD",
                                      "#5254A3", "#6B6ECF", "#9C9EDE", "#637939", "#B5CF6B",
                                      "#CEDB9C", "#8CA252", "#BD9E39", "#E7BA52", "#E7CB94",
                                      "#843C39", "#AD494A", "#D6616B", "#E7969C", "#7B4173",
                                      "#A55194", "#CE6DBD", "#DE9ED6", "#3182BD", "#6BAED6",
                                      "#9ECAE1", "#C6DBEF", "#E6550D", "#FD8D3C", "#FDAE6B",
                                      "#FDD0A2", "#31A354", "#74C476", "#A1D99B", "#C7E9C0",
                                      "#756BB1", "#9E9AC8", "#BCBDDC", "#DADAEB", "#636363",
                                      "#969696", "#BDBDBD", "#D9D9D9", "#393B79", "#5254A3",
                                      "#6B6ECF", "#9C9EDE", "#637939", "#8CA252", "#B5CF6B",
                                      "#CEDB9C", "#8C6D31", "#BD9E39", "#E7BA52", "#E7CB94"),
                           breaks = seq(0:24)-1) +
        scale_shape_manual(values = shape_levels) +
        scale_x_continuous(limits = c(minc_plot, maxc_plot), breaks = seq(minc_plot, maxc_plot, by = 2)) +
        scale_y_continuous(limits = c(minrt_plot, maxrt_plot)) +

        geom_point(
          data = subset(filter_data_plot, !is.na(comments) & comments != ""),
          aes(x = C, y = RTsec),
          shape = 22,
          size = 5,
          fill = "black",
          alpha = 0.25,
          color = NA,
          inherit.aes = FALSE
        )
    }

    # Add subtitle with lipid classes
    # Include which classes were used to build the models ......................
    models <- rtdata$filtering$models
    filter_info <- rtdata$filtering$information
    lipid_classes <- NULL
    if (!is.null(models[[lipidclass]]$model$data) && !is.null(models[[lipidclass]]$model$data$LipidSubclass)) {
      lipid_classes <- unique(models[[lipidclass]]$model$data$LipidSubclass)
    } else if (!is.null(models[[lipidclass]]$model$model) && !is.null(models[[lipidclass]]$model$model$LipidSubclass)) {
      lipid_classes <- unique(models[[lipidclass]]$model$model$LipidSubclass)
    } else {
    }

    subtitle_text <- character(0)
    if (!is.null(lipid_classes)) {
      if (length(lipid_classes) == 1) {
        subtitle_text <- paste0("This model was built using the lipid subclass: ", lipid_classes, ".")
      } else {
        subtitle_text <- paste0("This model was built using the following lipid subclasses: ",
                                paste0(lipid_classes, collapse = ", "), ".")
      }
    } else {
      subtitle_text <- "No model for this lipid subclass."
    }

    subtitle_text <- paste(
      subtitle_text, "\n ·",
      # kept = keep + recover
      (ifelse(is.na(filter_info$keep[filter_info$LipidSubclass == lipidclass]), 0,
              filter_info$keep[filter_info$LipidSubclass == lipidclass]) +
         ifelse(is.na(filter_info$recover[filter_info$LipidSubclass == lipidclass]), 0,
                filter_info$recover[filter_info$LipidSubclass == lipidclass])),
      "IDs kept,",
      # removed
      ifelse(is.na(filter_info$remove[filter_info$LipidSubclass == lipidclass]), 0,
             filter_info$remove[filter_info$LipidSubclass == lipidclass]),
      "IDs removed."
    )


    p <- p +
      ggtitle(subtitle_text) +
      theme(plot.title = element_text(size = 9, color = "grey41"))


    # If there is a prediction_layer, plot it ..................................
    if (exists("prediction_layer")) {
      if (length(prediction_layer) != 0) {
        p <- p +
          prediction_layer +
          theme(legend.position = "none")
      }
    }

    # Generate the interactive plot from the ggplot ............................
    plotly::layout(
      plotly::ggplotly(
        p,
        tooltip = c("text"),
        source = "manualFilterPlot_source"
      ),
      dragmode = "lasso"

    )

  }


  #============================================================================#
  # -- Render the initial Manual Filter plot
  #============================================================================#
  output$manualFilterPlot <- plotly::renderPlotly({
    rtdata <- rtdataReactive()
    if (!is.null(rtdata$filtering$data)) {
      rtdata$filtering$data <-
        create_other_isobars(rtdata$filtering$data,
                             typeanalysis = "filtering")
    }
    rtdataReactive(rtdata)

    req(input$lipid_class_selector)
    renderFilterPlot()
  })










  #============================================================================#
  # -- "Filter Combined Results" graphic interface.
  #============================================================================#
  update_filter_combined_ui <- function(session) {

    rtdata <- rtdataReactive()

    # Extract DB values for each Lipid Class
    db_choices <- DBextraction(rtdata)
    dbChoicesReactive(db_choices) # The DBs of each lipid class
    dbSelectionsReactive(db_choices) # The DBs selected in each lipid class by the user

    # UI: Render content of the filter results tab .............................
    output$filter_combined_ui <- renderUI({
      tagList(
        # Lipid class selection and graph download
        div(
          style = "margin-top: 0px; border-bottom: 1px solid #ddd; padding-bottom: 0px;",
          div(
            class = "no-width-change-tabs",
            fluidRow(
              column(width = 6,
                     tagList(
                       tags$div(
                         style = "display: flex; justify-content: flex-end; gap: 10px; margin-top: 0px; margin-bottom: 0px;",

                         tags$div(
                           style = "display: flex; align-items: center; justify-content: center; height: 100%;",

                           # Tutorial Image
                           actionButton(
                             inputId = "manual_filterCombined_info_btn",
                             label = tags$img(
                               src = "iconohorizontal_buttoninfo.png",
                               style = "height: auto; width: auto; margin-right: 10px; height: 50px;"
                             ),
                             style = "background-color: transparent; border: none; padding: 0;
        margin-right: 15px; cursor: pointer;"
                           ),

                           # Floating menu initially hidden
                           tags$div(
                             id = "manual_filterCombined_dropdown_menu",
                             style = " display: none; position: fixed; top: 50%; left: 50%;
        transform: translate(-50%, -50%); height: 90vh; width: 80vw;
        overflow-y: auto; background-color: rgba(255, 255, 255, 1);
        padding: 20px; font-size: 12px; text-align: right; border: none;
        border-radius: 5px; box-shadow: 0 4px 20px rgba(0,0,0,0.3); z-index: 1000;",

                             # Zoom controls
                             tags$div(
                               style = "position: sticky; top: 0; background-color: rgba(255,255,255,0.95);
          z-index: 10; padding-bottom: 5px; text-align: right;",
                               actionButton(
                                 "zoom_out_dropdown_filterCombined", "-",
                                 style = "width: 25px; height: 25px; padding: 0; font-size: 16px; vertical-align: middle;"
                               ),
                               span(
                                 id = "zoom_level_dropdown_filterCombined", "100%",
                                 style = "font-weight: bold; margin: 0 8px; vertical-align: middle;"
                               ),
                               actionButton(
                                 "zoom_in_dropdown_filterCombined", "+",
                                 style = "width: 25px; height: 25px; padding: 0; font-size: 16px; vertical-align: middle;"
                               )
                             ),

                             # Filtering tutorial image
                             tags$img(
                               id = "manualFilterCombined_dropdown",
                               src = "info_ManualFilter.png",
                               style = "width: 90%; height: auto; display: block; margin: 0 auto;"
                             ),

                             # JavaScript script for zoom
                             tags$script(HTML("
      var zoomDropdown_filterCombined = 100;
      var imgDropdown_filterCombined = document.getElementById('manualFilterCombined_dropdown');
      var zoomLevelSpanDropdown_filterCombined = document.getElementById('zoom_level_dropdown_filterCombined');

      document.getElementById('zoom_in_dropdown_filterCombined').onclick = function() {
        if (zoomDropdown_filterCombined < 200) {
          zoomDropdown_filterCombined += 10;
          imgDropdown_filterCombined.style.width = zoomDropdown_filterCombined + '%';
          zoomLevelSpanDropdown_filterCombined.textContent = zoomDropdown_filterCombined + '%';
        }
      };
      document.getElementById('zoom_out_dropdown_filterCombined').onclick = function() {
        if (zoomDropdown_filterCombined > 30) {
          zoomDropdown_filterCombined -= 10;
          imgDropdown_filterCombined.style.width = zoomDropdown_filterCombined + '%';
          zoomLevelSpanDropdown_filterCombined.textContent = zoomDropdown_filterCombined + '%';
        }
      };
    "))
                           ),

                           # Script to show/hide the dropdown when clicking the circle
                           tags$script(HTML("
          var manualFilterCombinedBtn = document.getElementById('manual_filterCombined_info_btn');
          var manualFilterCombinedMenu = document.getElementById('manual_filterCombined_dropdown_menu');

          manualFilterCombinedBtn.onclick = function(event) {
            event.stopPropagation();
            if (manualFilterCombinedMenu.style.display === 'none') {
              manualFilterCombinedMenu.style.display = 'block';
            } else {
              manualFilterCombinedMenu.style.display = 'none';
            }
          };

          document.addEventListener('click', function(event) {
            if (!manualFilterCombinedMenu.contains(event.target) && !manualFilterCombinedBtn.contains(event.target)) {
              manualFilterCombinedMenu.style.display = 'none';
            }
          });
        "))
                         ),

                       )
                     )
              )

            ),
            # Tabs for results and manual filtering
            tabsetPanel(
              type = "tabs",

              # Tab: Manual Filter
              tabPanel(
                "Manual Combined Filter",
                div(
                  style = "margin-top: 20px; margin-bottom: 40px;",
                  fluidRow(
                    column(2,

                           # Lipid Subclass selector
                           tags$details(
                             id = "subclasses_details_filter_combined",
                             style = "margin-bottom: 20px; border: 1px solid #dcdcdc;
                             border-radius: 6px; background-color: #fafafa; padding: 8px 12px;",

                             tags$summary(
                               style = "cursor: pointer; outline: none;",

                               tags$div(
                                 style = "font-size: 14px; font-weight: bold; color: #333; margin-bottom: 4px;",
                                 # icon("database"),
                                 # icon("flask"),
                                 icon("layer-group"),
                                 HTML("&nbsp;Lipid Subclass")
                               ),

                               tags$div(
                                 style = "font-size: 11px; color: #555; margin-top: 2px; text-align: center;",
                                 HTML("Choose which subclasses are displayed")
                               ),

                               tags$div(
                                 style = "text-align: center; font-size: 10px; color: #777; margin-top: 5px;",
                                 "click to expand / collapse"
                               )
                             ),

                             # SOLO el contenido interno es reactivo
                             div(
                               style = "margin-top: 15px;",
                               uiOutput("subclasses_filter_combined_content")
                             )
                           ),


                           tags$details(
                             style = "margin-bottom: 15px; padding: 8px 12px; border: 1px solid #dcdcdc;
                             border-radius: 6px; background-color: #fafafa;",

                             ## HEADER — todo clicable
                             tags$summary(
                               style = "cursor: pointer; outline: none;",

                               # Title
                               tags$div(
                                 style = "font-size: 14px; font-weight: bold; color: #333; margin-bottom: 4px;",
                                 icon("chart-area"),
                                 HTML("&nbsp;Graph style")
                               ),

                               # Explanation (always visible)
                               tags$div(
                                 style = "font-size: 11px; color: #555; margin-top: 2px; text-align: center; ",
                                 "Choose which data layers are displayed."
                               ),

                               # Click hint (below explanation)
                               tags$div(
                                 style = "text-align: center; font-size: 10px; color: #777; margin-top: 5px;",
                                 "click to expand / collapse"
                               )
                             ),

                             ## CONTENT — solo el select
                             # div(
                             #   style = "margin-top: 10px;",
                             #
                             #   selectInput(
                             #     inputId = "graph_style_filter_combined",
                             #     label   = NULL,   # quitamos el título interno
                             #     choices = list(
                             #       "Combined Data" = c(
                             #         "Kept + Predicted Data",
                             #         "Kept + Removed Data"
                             #       ),
                             #       "Single layers" = c(
                             #         "Kept Data",
                             #         "Predicted Data",
                             #         "Removed Data"
                             #       )
                             #     ),
                             #     selected = "Kept + Predicted Data",
                             #     width = "100%"
                             #   )
                             # )
                             div(
                               style = "margin-top: 10px;",

                               selectInput(inputId = "graph_style_filter_combined",
                                           label   = NULL,

                                           # Dynamically generate choices depending on prediction data availability
                                           choices = {

                                             # Check if prediction candidates exist and are not empty
                                             has_prediction <- (!is.null(rtdataReactive()$prediction$candidates) &&
                                                                  nrow(rtdataReactive()$prediction$candidates) > 0)

                                             # Choices when prediction data is available
                                             if (has_prediction) {

                                               list("Combined Data" = c("Kept Data + Accepted Candidates",
                                                                        "Kept Data + Removed Data",
                                                                        "Kept Data + Removed Data + Accepted Candidates"),
                                                    "Single layers" = c("Kept Data",
                                                                        "Accepted Candidates",
                                                                        "Removed Data"))

                                             } else {

                                               # Choices when prediction data is not available
                                               list("Combined Data" = c("Kept Data + Removed Data"),
                                                    "Single layers" = c("Kept Data", "Removed Data"))
                                             }
                                           },

                                           # Default selected option
                                           selected = if (!is.null(rtdataReactive()$prediction$candidates) &&
                                                          nrow(rtdataReactive()$prediction$candidates) > 0) {
                                             "Kept Data + Accepted Candidates"
                                           } else {
                                             # "Kept Data + Removed Data"
                                             "Kept Data"
                                           },

                                           width = "100%"
                               )
                             )
                           ),

                           # RT and C filtering (stacked, collapsible)
                           tags$details(
                             style = "margin-bottom: 15px; border: 1px solid #dcdcdc; border-radius: 6px;
                             background-color: #fafafa; padding: 8px 12px;",

                             ## Header
                             tags$summary(
                               style = "display: block; cursor: pointer; outline: none;",

                               # Línea superior: icono + título
                               tags$div(
                                 style = "display: flex; align-items: center; font-size: 14px;
                                 font-weight: bold; color: #333; margin-bottom: 2px;",
                                 icon("filter"),
                                 HTML("&nbsp;&nbsp;RT - C Ranges")
                               ),

                               # Texto explicativo (siempre visible)
                               tags$div(
                                 style = "font-size: 11px; color: #555; line-height: 1.3; margin-bottom: 2px; text-align: center;",
                                 HTML(
                                   "Use these controls to modify the Retention Time and Carbon Number ranges displayed."
                                 )
                               ),

                               # Click hint CENTRADO debajo
                               tags$div(
                                 style = "text-align: center; font-size: 10px; color: #777; margin-top: 5px;",
                                 "click to expand / collapse"
                               )
                             ),

                             ## Content (STACKED)
                             div(
                               style = "margin-top: 10px;",

                               # RT range
                               div(
                                 style = "padding: 4px; border: 1px dashed #d0d0d0; border-radius: 4px;
                                 background-color: #ffffff; margin-bottom: 10px;",
                                 uiOutput("rt_range_filter_combined_ui")
                               ),

                               # C range
                               div(
                                 style = "padding: 4px; border: 1px dashed #d0d0d0; border-radius: 4px;
                                 background-color: #ffffff;",
                                 uiOutput("c_range_filter_combined_ui")
                               )
                             )
                           ),

                           # DB selector
                           tags$details(
                             id = "db_details_filter_combined",
                             style = "margin-bottom: 20px; border: 1px solid #dcdcdc;
                             border-radius: 6px; background-color: #fafafa; padding: 8px 12px;",

                             tags$summary(
                               style = "cursor: pointer; outline: none;",

                               tags$div(
                                 style = "font-size: 14px; font-weight: bold; color: #333; margin-bottom: 4px;",
                                 icon("database"),
                                 HTML("&nbsp;DB values")
                               ),

                               tags$div(
                                 style = "font-size: 11px; color: #555; margin-top: 2px; text-align: center;",
                                 HTML("Select the DB values affected by <b>Keep</b> and <b>Remove</b> actions.")
                               ),

                               tags$div(
                                 style = "text-align: center; font-size: 10px; color: #777; margin-top: 5px;",
                                 "click to expand / collapse"
                               )
                             ),

                             # SOLO el contenido interno es reactivo
                             div(
                               style = "margin-top: 15px;",
                               uiOutput("db_filter_combined_content")
                             )
                           ),
                           div(
                             style = "margin-top: 15px; padding: 10px; border: 1px solid #dcdcdc;
                             border-radius: 6px; background-color: #fafafa;",

                             # Title
                             tags$div(
                               style = "font-size: 14px; font-weight: bold; margin-bottom: 4px;",
                               icon("sliders-h"),
                               HTML("&nbsp;Modify filtering and predictions")
                             ),

                             # Keep button
                             actionButton(
                               "keepFilterButton_combined",
                               label = HTML('<span style="font-size:16px; margin-right:6px;">&#x25B2;</span><b style="font-size:14px;">Keep</b>'),
                               style = "width: 100%; margin-top: 6px; margin-bottom: 6px; background-color: #ecf0f1;
                               color: black; border: 1px solid #d0d0d0; padding: 4px 8px;"),

                             # Remove button
                             actionButton(
                               "removeFilterButton_combined",
                               label = HTML('<span style=\"font-size:16px; margin-right:6px;\">&#x25B3;</span><b style=\"font-size:14px;\">Remove</b>'),
                               style = "width: 100%; margin-bottom: 8px; background-color: #ecf0f1;
                               color: black; border: 1px solid #d0d0d0; padding: 4px 8px;"
                             ),

                             uiOutput("saveFilterButtonUI_combined"),

                             # Save hint text
                             tags$div(
                               style = "font-size: 11px; color: #666; margin-top: 2px; text-align: center;",
                               icon("info-circle"),
                               HTML("&nbsp;Changes are only applied after saving")
                             )
                           )
                    ),

                    column(8,

                           # Interactive Prediction Graph
                           div(
                             class = "full-height",
                             plotly::plotlyOutput("manualFilterPlot_combined", height = "700px")
                           ),

                           ## Info box under the graph
                        #    div(
                        #      style = "margin-top: 12px; padding: 10px 12px; border: 1px solid #dcdcdc;
                        # border-left: 4px solid #6c757d; border-radius: 4px; background-color: #f8f9fa;
                        # font-size: 12px; color: #333;",
                        #      icon("info-circle"),
                        #      tags$strong(" Filtering impact: "),
                        #      "If you modify the filtering, the predictions will be recalculated accordingly."
                        #    )
                    ),
                    column(2,
                           # Legend image of the filter graph. If you click on the image, it doubles in size. To return it to its original size, tap anywhere on the screen.
                           div(
                             class = "full-width",
                             tags$img(id = "zoomable-image_filter_combined",src = "legend_manualFilter_combined.png",
                                      width = "90%",
                                      style = "cursor: pointer; transition: transform 0.3s ease-in-out;
                          position: relative; top: 0; left: 0;"
                             ),
                             tags$script(HTML("
                            $(document).ready(function() {
                              $('#zoomable-image_filter_combined').click(function(event) {
                                event.stopPropagation();
                                $(this).toggleClass('zoomed');
                              });
                              $(document).click(function() {
                                $('#zoomable-image_filter_combined').removeClass('zoomed');
                              });
                            });
                          ")),
                             tags$style(HTML("
                            #zoomable-image_filter_combined {
                              max-width: 100%;
                            }
                            #zoomable-image_filter_combined.zoomed {
                              transform: scale(2);
                              transform-origin: top left;
                              position: absolute;
                              top: 0;
                              left: 0;
                            }
                        "))
                           ),

                           div(
                             style = "margin-top: 20px; padding: 10px; border: 1px solid #dcdcdc;
                             border-radius: 6px; background-color: #fafafa;",

                             # Title
                             tags$div(
                               style = "font-size: 14px; font-weight: bold; margin-bottom: 4px;",
                               icon("highlighter"),
                               HTML("&nbsp;Highlight data")
                             ),

                             # Subtitle / explanation
                             tags$div(
                               style = "font-size: 11px; color: #555; margin-bottom: 8px; text-align: center;",
                               "Mark or Unmark lipids in the plot to highlight them."
                             ),

                             # Mark button
                             actionButton(
                               "markFilterButton_combined",
                               label = HTML(paste(icon("tag"), '<b style="font-size: 13px; margin-left: 6px;">Mark</b>')),
                               style = "width: 100%; margin-bottom: 6px; background-color: #ecf0f1;
                               color: black; border: 1px solid #d0d0d0; padding: 4px 8px;"),

                             # Unmark button
                             actionButton(
                               "unmarkFilterButton_combined",
                               label = HTML(paste(icon("times"),'<b style="font-size: 13px; margin-left: 6px;">Unmark</b>')),
                               style = "width: 100%; background-color: #ecf0f1; color: black;
                               border: 1px solid #d0d0d0; padding: 4px 8px; margin-bottom: 4px;"),

                             # Save mark button
                             actionButton(
                               "saveMarkFilterButton_combined",
                               label = HTML(paste(icon("save"),'<b style="font-size: 13px; margin-left: 6px;">Save</b>')),
                               style = "width: 100%; background-color: #ecf0f1; color: black;
                               border: 1px solid #d0d0d0; padding: 4px 8px; margin-bottom: 4px;"),

                             # Save hint text
                             tags$div(
                               style = "font-size: 11px; color: #666; margin-top: 2px; text-align: center;",
                               icon("info-circle"),
                               HTML("&nbsp;Changes are only applied after saving")
                             )
                           )

                    )

                    )
                )
              )

            )
          )
        )
      )
    })

    output$saveFilterButtonUI_combined <- renderUI({

      pending_changes <- (!is.null(selectedindices_filter_combined()) &&
                            length(selectedindices_filter_combined()) > 0) ||
        (!is.null(selectedindices_prediction_combined()) &&
           length(selectedindices_prediction_combined()) > 0)

      if (pending_changes) {

        actionButton(
          "saveFilterButton_combined",
          label = HTML(
            paste(
              icon("save"),
              '<b style="font-size:14px; margin-left:6px;">Save</b>'
            )
          ),
          style = "
        width: 100%;
        margin-bottom: 4px;
        background-color: #d6ebff;
        color: #003b6f;
        border: 1px solid #9bc7ef;
        padding: 4px 8px;
        font-weight: bold;
      "
        )

      } else {

        actionButton(
          "saveFilterButton_combined",
          label = HTML(
            paste(
              icon("save"),
              '<b style="font-size:14px; margin-left:6px;">Save</b>'
            )
          ),
          style = "
        width: 100%;
        margin-bottom: 4px;
        background-color: #ecf0f1;
        color: black;
        border: 1px solid #d0d0d0;
        padding: 4px 8px;
      "
        )

      }
    })



    # Filtered Lipid Subclasses Interface Values ...............................
    observeEvent(lipid_classes_combined(), {

      current_classes <- lipid_classes_combined()
      req(current_classes)

      old_selection <- selected_lipid_classes_combined()

      # Mantener las que el usuario ya tenía seleccionadas
      still_valid <- intersect(old_selection, current_classes)

      # Añadir solo clases nuevas
      new_classes <- setdiff(current_classes, old_selection)

      updated_selection <- c(still_valid, new_classes)

      selected_lipid_classes_combined(updated_selection)

    }, ignoreInit = FALSE)


    # Render de checkboxes -----------------------------------------------------

    output$subclasses_filter_combined_content <- renderUI({

      current_classes <- lipid_classes_combined()
      req(current_classes)

      checkboxGroupInput(
        inputId = "subclasses_values_filter_combined",
        label = NULL,
        choices = current_classes,
        selected = selected_lipid_classes_combined()
      )

    })


    # Guardar cambios del usuario ----------------------------------------------

    observeEvent(input$subclasses_values_filter_combined, {

      selected_lipid_classes_combined(
        input$subclasses_values_filter_combined
      )

    }, ignoreInit = TRUE)


    # # Filtered DB Interface Values ---------------------------------------------
    #
    # output$db_filter_combined_content <- renderUI({
    #
    #   req(dbChoicesReactive())
    #   req(dbSelectionsReactive())
    #   req(selected_lipid_classes_combined())
    #
    #   selected_classes <- selected_lipid_classes_combined()
    #
    #   if (length(selected_classes) == 0) {
    #     return(NULL)
    #   }
    #
    #
    #   # Obtener todos los DB values de todas las clases seleccionadas
    #   choices_vec <- unlist(
    #     dbChoicesReactive()$filtering[selected_classes]
    #   )
    #
    #   choices_vec <- as.numeric(choices_vec)
    #   choices_vec <- sort(unique(choices_vec[!is.na(choices_vec)]))
    #
    #
    #   # Obtener todas las selecciones previas de esas clases
    #   previous_selected <- unlist(
    #     dbSelectionsReactive()$filtering[selected_classes]
    #   )
    #
    #   selected_db_values <- intersect(
    #     unique(previous_selected),
    #     choices_vec
    #   )
    #
    #
    #   div(
    #     checkboxGroupInput(
    #       "db_values_filter_combined",
    #       label = NULL,
    #       choices = choices_vec,
    #       selected = selected_db_values,
    #       inline = FALSE
    #     ),
    #
    #     if (length(choices_vec) == 0) {
    #       tags$div(
    #         style = "font-size: 11px; color: #777; margin-top: 6px;",
    #         "No DB values available"
    #       )
    #     }
    #   )
    #
    # })
    #
    #
    # # Save selection DB for all active classes .................................
    # observeEvent(input$db_values_filter_combined, {
    #   req(selected_lipid_classes_combined())
    #
    #   selected_classes <- selected_lipid_classes_combined()
    #
    #   db_selections <- dbSelectionsReactive()
    #
    #   # Aplicar la misma selección a todas las clases activas
    #   for (lipid_class in selected_classes) {
    #     db_selections$filtering[[lipid_class]] <-input$db_values_filter_combined
    #   }
    #
    #   dbSelectionsReactive(db_selections)
    # })

    # Filtered DB Interface Values ---------------------------------------------

    output$db_filter_combined_content <- renderUI({

      req(dbChoicesReactive())
      req(dbSelectionsReactive())
      req(selected_lipid_classes_combined())

      selected_classes <- selected_lipid_classes_combined()

      if (length(selected_classes) == 0) {
        return(NULL)
      }

      # Get all DB values from selected lipid classes
      choices_vec <- unlist(
        dbChoicesReactive()$filtering[selected_classes]
      )

      # If graph style contains "Predicted",
      # also include prediction DB values
      if (!is.null(input$graph_style_filter_combined) &&
        grepl("Candidates", input$graph_style_filter_combined)) {

        prediction_vec <- unlist(
          dbChoicesReactive()$prediction[selected_classes]
        )

        choices_vec <- c(choices_vec, prediction_vec)
      }

      # Clean and sort DB values
      choices_vec <- as.numeric(choices_vec)
      choices_vec <- sort(unique(choices_vec[!is.na(choices_vec)]))

      # Get previous selections from selected classes
      previous_selected <- unlist(
        dbSelectionsReactive()$filtering[selected_classes]
      )

      # Keep only valid selected values
      selected_db_values <- intersect(
        unique(previous_selected),
        choices_vec
      )

      div(
        checkboxGroupInput(
          "db_values_filter_combined",
          label = NULL,
          choices = choices_vec,
          selected = selected_db_values,
          inline = FALSE
        ),

        # Message when no DB values are available
        if (length(choices_vec) == 0) {
          tags$div(
            style = "font-size: 11px; color: #777; margin-top: 6px;",
            "No DB values available"
          )
        }
      )

    })


    # Save DB selections for all active classes -------------------------------

    observeEvent(input$db_values_filter_combined, {
      req(selected_lipid_classes_combined())

      selected_classes <- selected_lipid_classes_combined()
      db_selections <- dbSelectionsReactive()

      # Apply the same selection to all active lipid classes
      for (lipid_class in selected_classes) {
        db_selections$filtering[[lipid_class]] <- input$db_values_filter_combined
      }

      dbSelectionsReactive(db_selections)
    })


    # Filtered C Range Interface Values ........................................
    output$c_range_filter_combined_ui <- renderUI({
      req(minC_filter_combined())
      req(maxC_filter_combined())
      c_min <- minC_filter_combined()
      c_max <- maxC_filter_combined()

      sliderInput("c_range_filter_combined", "CN Range:", min = c_min, max = c_max,
                  value = c(c_min, c_max), step = 1, width = "80%")
    })


    # Filtered RT Range Interface Values .......................................
    output$rt_range_filter_combined_ui <- renderUI({
      req(minRT_filter_combined(), maxRT_filter_combined())

      rt_min <- minRT_filter_combined()
      rt_max <- maxRT_filter_combined()

      div(
        style = "margin-bottom: 0px;",

        tags$div(
          style = "font-size: 12px; font-weight: bold; margin-bottom: 2px; text-align: center;",
          "Retention Time Range"
        ),

        sliderInput(
          "rt_range_filter_combined",
          label = NULL,
          min = rt_min,
          max = rt_max,
          value = c(rt_min, rt_max),
          step = 10,
          width = "100%"
        )
      )
    })


    # Filtered C Range Interface Values ......................................
    output$c_range_filter_combined_ui <- renderUI({
      req(minC_filter_combined())
      req(maxC_filter_combined())
      c_min <- minC_filter_combined()
      c_max <- maxC_filter_combined()

      div(
        style = "margin-bottom: 0px;",

        tagList(
          tags$div(
            style = "font-size: 12px; font-weight: bold; margin-bottom: 2px; text-align: center;",
            "Carbon Number Range"
          ),
          sliderInput(
            "c_range_filter_combined",
            label = NULL,
            min = c_min,
            max = c_max,
            value = c(c_min, c_max),
            step = 1,
            width = "100%"
          )))
    })
  }




  #============================================================================#
  # -- 'Manual Combiend' Graph Controls
  #============================================================================#
  # REMOVE COMBINED ............................................................
  observeEvent(input$removeFilterButton_combined, {

    # Modified lipid indices in the interactive filtering graph
    selected <- plotly::event_data("plotly_selected", source = "manualFilterPlot_combined_source")

    rtdata <- rtdataReactive()

    # Keep the indices of selected lipids that match the user-selected DBs
    if (!is.null(selected)) {

      new_indices <- as.numeric(unique(selected$customdata))

      if (!is.null(input$graph_style_filter_combined) &&
          grepl("Data", input$graph_style_filter_combined) &&
          !is.null(rtdata$filtering$data) &&
          nrow(rtdata$filtering$data) > 0) {

        indices_filter <- intersect(
          new_indices,
          rtdata$filtering$data$index[
            rtdata$filtering$data$DB %in% input$db_values_filter_combined
          ]
        )

        if (!is.null(indices_filter) & length(indices_filter) > 0) {

          if ("Filter_final" %in% colnames(rtdata$filtering$data)) {

            rtdata$filtering$data$Filter_final[
              rtdata$filtering$data$index %in% indices_filter
            ] <- "Removed manually by the user"

          } else {

            idx <- rtdata$filtering$data$index %in% indices_filter

            rtdata$filtering$data$checkfilter[idx] <- ifelse(

              ## PRIORITY 1: Removed with RT modelling (keep it as is)
              rtdata$filtering$data$checkmodel[idx] == "Removed with RT modelling",

              "Removed with RT modelling",

              ifelse(

                ## PRIORITY 2: manual removal cases
                rtdata$filtering$data$Filter[idx] %in% c("Kept") |

                  (rtdata$filtering$data$checkmodel[idx] %in% c("Recovered with RT modelling") &
                     rtdata$filtering$data$Filter[idx] %in% c(
                       "Confidence level among isobars",
                       "Removed by RT rules",
                       "In-source fragments",
                       "Removed distant data",
                       "Outlier",
                       "RT coherence between subclasses")),

                "Removed manually by the user",

                ## DEFAULT
                rtdata$filtering$data$Filter[idx]
              )
            )
          }

          selected_indices <- selectedindices_filter_combined()
          selected_indices <- c(selected_indices, indices_filter)
          selectedindices_filter_combined(selected_indices)

          # Render the updated interactive graph
          output$manualFilterPlot <- plotly::renderPlotly({
            req(input$lipid_class_selector)
            renderFilterPlot()
          })

        }
      }



      if (!is.null(input$graph_style_filter_combined) &&
          grepl("Candidates", input$graph_style_filter_combined) &&
          !is.null(rtdata$prediction$candidates) &&
          nrow(rtdata$prediction$candidates) > 0) {

        indices_prediction <- intersect(
          new_indices,
          rtdata$prediction$candidates$index[
            rtdata$prediction$candidates$DB %in% input$db_values_filter_combined
          ]
        )

        if (!is.null(indices_prediction) & length(indices_prediction) > 0) {

          if ("Filter_final" %in% colnames(rtdata$prediction$candidates)) {

            # Manually assign 'Removed manually by the user' to the checkcandidate field
            rtdata$prediction$candidates$Filter_final[
              rtdata$prediction$candidates$index %in% indices_prediction
            ] <- "Removed manually by the user"

          } else {
            # Manually assign 'Removed manually by the user' to the checkcandidate field
            rtdata$prediction$candidates$checkcandidate[
              rtdata$prediction$candidates$index %in% indices_prediction
            ] <-ifelse(
              rtdata$prediction$candidates$checkscore_rtpred[
                rtdata$prediction$candidates$index %in% indices_prediction
              ] %in% c("Kept"),
              "Removed manually by the user",
              rtdata$prediction$candidates$checkscore_rtpred[
                rtdata$prediction$candidates$index %in% indices_prediction]
            )
          }

          selected_indices <- selectedindices_prediction_combined()
          selected_indices <- c(selected_indices, indices_prediction)
          selectedindices_prediction_combined(selected_indices)

          # Render the updated interactive graph
          output$manualPredictionPlot <- plotly::renderPlotly({
            renderPredictionGraph()
          })

        }
      }

      # Update rtdataReactive
      rtdataReactive(rtdata)


      output$manualFilterPlot_combined <- plotly::renderPlotly({

        rtdata <- rtdataReactive()

        if (!is.null(rtdata$filtering$data)) {
          rtdata$filtering$data <- create_other_isobars(rtdata$filtering$data,
                                                        typeanalysis = "filtering")
        }
        if (!is.null(rtdata$prediction$candidates)) {
          rtdata$prediction$candidates <- create_other_isobars(rtdata$prediction$candidates,
                                                               typeanalysis = "prediction")
        }

        rtdataReactive(rtdata)

        req(lipid_classes_combined())
        req(selected_lipid_classes_combined())

        lipid_classes_combined()
        selected_lipid_classes_combined()

        renderFilterPlot_combined()
      })

    }
  })

  # KEEP COMBINED ..............................................................
  observeEvent(input$keepFilterButton_combined, {

    # Modified lipid indices in the interactive filtering graph
    selected <- plotly::event_data("plotly_selected", source = "manualFilterPlot_combined_source")

    rtdata <- rtdataReactive()

    # Keep the indices of selected lipids that match the user-selected DBs
    if (!is.null(selected)) {

      new_indices <- as.numeric(unique(selected$customdata))

      if (!is.null(input$graph_style_filter_combined) &&
          grepl("Data", input$graph_style_filter_combined) &&
          !is.null(rtdata$filtering$data) &&
          nrow(rtdata$filtering$data) > 0) {

        indices_filter <- intersect(
          new_indices,
          rtdata$filtering$data$index[
            rtdata$filtering$data$DB %in% input$db_values_filter_combined
          ]
        )

        if (!is.null(indices_filter) & length(indices_filter) > 0) {

          if ("Filter_final" %in% colnames(rtdata$filtering$data)) {

            rtdata$filtering$data$Filter_final[
              rtdata$filtering$data$index %in% indices_filter
            ] <- "Kept manually by the user"

          } else {

            idx <- rtdata$filtering$data$index %in% indices_filter

            checkmodel <- rtdata$filtering$data$checkmodel[idx]
            Filter     <- rtdata$filtering$data$Filter[idx]

            rtdata$filtering$data$checkfilter[idx] <- ifelse(

              ## PRIORITY 1: keep "Recovered with RT modelling"
              checkmodel == "Recovered with RT modelling",
              "Recovered with RT modelling",

              ifelse(

                ## PRIORITY 2: manual keep cases
                Filter %in% c("Confidence level among isobars", "Removed by RT rules",
                              "In-source fragments", "Removed distant data",
                              "Outlier", "RT coherence between subclasses") |
                  (checkmodel == "Removed with RT modelling" & Filter == "Kept"),

                "Kept manually by the user",

                ## DEFAULT
                Filter
              )
            )
          }

          selected_indices <- selectedindices_filter_combined()
          selected_indices <- c(selected_indices, indices_filter)
          selectedindices_filter_combined(selected_indices)

          # Render the updated interactive graph
          output$manualFilterPlot <- plotly::renderPlotly({
            req(input$lipid_class_selector)
            renderFilterPlot()
          })

        }
      }



      if (!is.null(input$graph_style_filter_combined) &&
          grepl("Candidates", input$graph_style_filter_combined) &&
          !is.null(rtdata$prediction$candidates) &&
          nrow(rtdata$prediction$candidates) > 0) {

        indices_prediction <- intersect(
          new_indices,
          rtdata$prediction$candidates$index[
            rtdata$prediction$candidates$DB %in% input$db_values_filter_combined
          ]
        )

        if (!is.null(indices_prediction) & length(indices_prediction) > 0) {

          if ("Filter_final" %in% colnames(rtdata$prediction$candidates)) {

            # Manually assign 'Removed manually by the user' to the checkcandidate field
            rtdata$prediction$candidates$Filter_final[
              rtdata$prediction$candidates$index %in% indices_prediction
            ] <- "Removed manually by the user"

          } else {
            # Manually assign 'Removed manually by the user' to the checkcandidate field
            rtdata$prediction$candidates$checkcandidate[
              rtdata$prediction$candidates$index %in% indices_prediction
            ] <-ifelse(
              rtdata$prediction$candidates$checkscore_rtpred[
                rtdata$prediction$candidates$index %in% indices_prediction
              ] %in% c("Kept"),
              "Removed manually by the user",
              rtdata$prediction$candidates$checkscore_rtpred[
                rtdata$prediction$candidates$index %in% indices_prediction]
            )
          }

          selected_indices <- selectedindices_prediction_combined()
          selected_indices <- c(selected_indices, indices_prediction)
          selectedindices_prediction_combined(selected_indices)

          # Render the updated interactive graph
          output$manualPredictionPlot <- plotly::renderPlotly({
            renderPredictionGraph()
          })

        }
      }

      # Update rtdataReactive
      rtdataReactive(rtdata)


      output$manualFilterPlot_combined <- plotly::renderPlotly({

        rtdata <- rtdataReactive()

        if (!is.null(rtdata$filtering$data)) {
          rtdata$filtering$data <- create_other_isobars(rtdata$filtering$data,
                                                        typeanalysis = "filtering")
        }
        if (!is.null(rtdata$prediction$candidates)) {
          rtdata$prediction$candidates <- create_other_isobars(rtdata$prediction$candidates,
                                                               typeanalysis = "prediction")
        }

        rtdataReactive(rtdata)

        req(lipid_classes_combined())
        req(selected_lipid_classes_combined())

        lipid_classes_combined()
        selected_lipid_classes_combined()

        renderFilterPlot_combined()
      })

    }
  })


  # SAVE FILTER ................................................................
  observeEvent(input$saveFilterButton_combined, {
    notif_id <- showNotification(HTML("The changes from the combined plot are being saved...  <span class='loader'></span>"),
                                 type = "default", duration = NULL)

    # Save the changes in the reactive rtdata
    rtdata <- rtdataReactive()

    indices_filter <- selectedindices_filter_combined()
    indices_prediction <- selectedindices_prediction_combined()

    # Changes in filtering .....................................................
    if (!is.null(indices_filter) & length(indices_filter) > 0) {
    # Rebuild models after saving the changes in filtering. (si no se ha hecho las predicciones aún)

      if (!"Filter_final" %in% colnames(rtdata$filtering$data)) {
        rtdata <- LipidRT::checkremovedids(rtdata = rtdata,
                                           distribution = input$distribution_filter,
                                           rtscorethr = input$score_filter)
      }

      rtdata <- LipidRT::organizertdata(rtdata = rtdata,
                                        typeanalysis = "filtering")

      # Redo the graphs with the changes
      rtdata <- LipidRT::filteringgraph(rtdata = rtdata)

      rtdata$filtering$data <- create_other_isobars(rtdata$filtering$data,
                                                    typeanalysis = "filtering")

      # Render automatic graphs
      output$filterPlot <- renderPlot({
        req(input$lipid_class_selector)
        lipid_class <- input$lipid_class_selector

        rtdata$filtering$plots[[lipid_class]]
      })

      # Render the updated interactive graph
      output$manualFilterPlot <- plotly::renderPlotly({
        req(input$lipid_class_selector)

        renderFilterPlot()
      })


      if (!is.null(rtdata$prediction$candidates)) {
        rtdata <- LipidRT::organizertdata(rtdata = rtdata,
                                          typeanalysis = "prediction")

        # Redo the graphs with the changes
        rtdata <- LipidRT::predictgraph(rtdata = rtdata,
                                        eval_predict_plot = TRUE)

        rtdata$prediction$candidates <- create_other_isobars(rtdata$prediction$candidates,
                                                             typeanalysis = "prediction")

        # Render automatic graphs
        output$predictionPlot <- renderPlot({
          req(input$lipid_class_selector)
          lipid_class <- input$lipid_class_selector

          rtdata$prediction$plots[[lipid_class]]
        })

        # Render the updated interactive graph
        output$manualPredictionPlot <- plotly::renderPlotly({
          req(input$lipid_class_selector)

          renderPredictionGraph()
        })

      }

    }


    if (!is.null(indices_prediction) & length(indices_prediction) > 0) {

      rtdata <- predictgraph(rtdata = rtdata,
                             classes = NULL,
                             eval_predict_plot = TRUE)

      rtdata <- organizertdata(rtdata = rtdata,
                               typeanalysis = "prediction")

      rtdata$filtering$data <- create_other_isobars(rtdata$filtering$data,
                                                    typeanalysis = "filtering")

      rtdata$prediction$candidates <- create_other_isobars(rtdata$prediction$candidates,
                                                           typeanalysis = "prediction")

      # Render automatic graphs
      output$predictionPlot <- renderPlot({
        req(input$lipid_class_selector)
        lipid_class <- input$lipid_class_selector

        rtdata$prediction$plots[[lipid_class]]
      })

      # Render the updated interactive graph
      output$manualPredictionPlot <- plotly::renderPlotly({
        req(input$lipid_class_selector)

        renderPredictionGraph()
      })

    }

    output$manualFilterPlot_combined <- plotly::renderPlotly({

      rtdata <- rtdataReactive()

      if (!is.null(rtdata$filtering$data)) {
        rtdata$filtering$data <- create_other_isobars(rtdata$filtering$data,
                                                      typeanalysis = "filtering")
      }
      if (!is.null(rtdata$prediction$candidates)) {
        rtdata$prediction$candidates <- create_other_isobars(rtdata$prediction$candidates,
                                                             typeanalysis = "prediction")
      }

      rtdataReactive(rtdata)

      req(lipid_classes_combined())
      req(selected_lipid_classes_combined())

      lipid_classes_combined()
      selected_lipid_classes_combined()


      renderFilterPlot_combined()
    })

    selectedindices_filter_combined(NULL)
    selectedindices_prediction_combined(NULL)

    rtdataReactive(rtdata)

    # Update feature table
    update_features_ui(session)

    # Notification: Saved successfully (lasts 3 seconds) .......................
    removeNotification(notif_id)  # Remove the 'saving' notification
    showNotification("The changes made in the combined plot have been successfully saved.",
                     type = "message", duration = 3)

  })





  # # MARK FILTER COMBINED .....................................................
  observeEvent(input$markFilterButton_combined, {

    showModal(modalDialog(
      title = NULL,

      div(
        style = "
        padding: 15px;
        background-color: #f8f9fa;
        border-radius: 8px;
      ",

        h4(icon("tag"), " Mark selected lipids"),
        p(
          style = "color: #666;",
          "Enter a comment that will be attached to all selected lipids."
        ),

        textAreaInput(
          "mark_text_combined",
          label = NULL,
          placeholder = "Write your comment here...",
          rows = 3
        )
      ),

      footer = tagList(
        modalButton("Cancel"),
        actionButton(
          "confirm_mark_filter_combined",
          "Apply comment",
          class = "btn-primary"
        )
      ),

      easyClose = TRUE
    ))
  })

  observe({
    shinyjs::toggleState(
      "confirm_mark_filter_combined",
      condition = nzchar(trimws(input$mark_text_combined))
    )
  })


  observeEvent(input$confirm_mark_filter_combined, {

    removeModal()

    selected <- plotly::event_data("plotly_selected", source = "manualFilterPlot_combined_source")

    rtdata <- rtdataReactive()

    if (!is.null(selected)) {

      new_indices <- as.numeric(unique(selected$customdata))

        indices_filter <- intersect(
          new_indices,
          rtdata$filtering$data$index[
            rtdata$filtering$data$DB %in% input$db_values_filter_combined
          ]
        )


      if (!is.null(indices_filter) & length(indices_filter) > 0) {
        # Usar el texto escrito por el usuario
        rtdata$filtering$data$comments[rtdata$filtering$data$index %in% indices_filter] <- input$mark_text_combined
      }


        if (!is.null(input$graph_style_filter_combined) &&
            grepl("Candidates", input$graph_style_filter_combined) &&
            !is.null(rtdata$prediction$candidates) &&
            nrow(rtdata$prediction$candidates) > 0) {

        indices_prediction <- intersect(
          new_indices,
          rtdata$prediction$candidates$index[
            rtdata$prediction$candidates$DB %in% input$db_values_filter_combined
          ]
        )

        if (!is.null(indices_prediction) & length(indices_prediction) > 0) {
          # Usar el texto escrito por el usuario
          rtdata$prediction$candidates$comments[rtdata$prediction$candidates$index %in% indices_prediction] <- input$mark_text_combined
        }

      }

      rtdataReactive(rtdata)


      output$manualFilterPlot_combined <- plotly::renderPlotly({

        rtdata <- rtdataReactive()

        if (!is.null(rtdata$filtering$data)) {
          rtdata$filtering$data <- create_other_isobars(rtdata$filtering$data,
                                                        typeanalysis = "filtering")
        }
        if (!is.null(rtdata$prediction$candidates)) {
          rtdata$prediction$candidates <- create_other_isobars(rtdata$prediction$candidates,
                                                               typeanalysis = "prediction")
        }

        rtdataReactive(rtdata)

        req(lipid_classes_combined())
        req(selected_lipid_classes_combined())

        lipid_classes_combined()
        selected_lipid_classes_combined()

        renderFilterPlot_combined()
      })


    }
  })


  # UNMARK FILTER COMBINED .....................................................
  observeEvent(input$unmarkFilterButton_combined, {

    # Modified lipid indices in the interactive prediction graph
    selected <- plotly::event_data("plotly_selected", source = "manualFilterPlot_combined_source")

    rtdata <- rtdataReactive()

    # Keep the indices of selected lipids that match the user-selected DBs
    if (!is.null(selected)) {
      new_indices <- as.numeric(unique(selected$customdata))

        indices_filter <- intersect(
          new_indices,
          rtdata$filtering$data$index[
            rtdata$filtering$data$DB %in% input$db_values_filter_combined
          ]
        )

        if (!is.null(indices_filter) & length(indices_filter) > 0) {
          rtdata$filtering$data$comments[rtdata$filtering$data$index %in% indices_filter] <- ""
        }


        if (!is.null(input$graph_style_filter_combined) &&
            grepl("Candidates", input$graph_style_filter_combined) &&
            !is.null(rtdata$prediction$candidates) &&
            nrow(rtdata$prediction$candidates) > 0) {

          indices_prediction <- intersect(
            new_indices,
            rtdata$prediction$candidates$index[
              rtdata$prediction$candidates$DB %in% input$db_values_filter_combined
            ]
          )

          if (!is.null(indices_prediction) & length(indices_prediction) > 0) {
            # Usar el texto escrito por el usuario
            rtdata$prediction$candidates$comments[rtdata$prediction$candidates$index %in% indices_prediction] <- ""
          }
        }

        # Update 'rtdataReactive'
        rtdataReactive(rtdata)

        # Render the updated interactive graph
        output$manualFilterPlot_combined <- plotly::renderPlotly({

          rtdata <- rtdataReactive()

          if (!is.null(rtdata$filtering$data)) {
            rtdata$filtering$data <- create_other_isobars(rtdata$filtering$data,
                                                          typeanalysis = "filtering")
          }
          if (!is.null(rtdata$prediction$candidates)) {
            rtdata$prediction$candidates <- create_other_isobars(rtdata$prediction$candidates,
                                                                 typeanalysis = "prediction")
          }

          rtdataReactive(rtdata)

          req(lipid_classes_combined())
          req(selected_lipid_classes_combined())

          lipid_classes_combined()
          selected_lipid_classes_combined()

          renderFilterPlot_combined()
        })
    }
  })


  # SAVE MARK FILTER COMBINED ..................................................
  observeEvent(input$saveMarkFilterButton_combined, {
    notif_id <- showNotification(HTML("The changes from the Highlight are being saved...  <span class='loader'></span>"),
                                 type = "default", duration = NULL)

    # Save the changes in the reactive rtdata
    rtdata <- rtdataReactive()

    # Changes in filtering .....................................................
    rtdata <- LipidRT::organizertdata(rtdata = rtdata,
                                      typeanalysis = "filtering")

    # Redo the graphs with the changes
    rtdata <- LipidRT::filteringgraph(rtdata = rtdata)

    if (!is.null(rtdata$prediction$candidates)) {
      # Modify the feature matrix with the changes
      rtdata <- LipidRT::organizertdata(rtdata = rtdata,
                                        typeanalysis = "prediction")

      # Redo the graphs with the changes
      rtdata <- LipidRT::predictgraph(rtdata = rtdata,
                                      eval_predict_plot = TRUE)
    }

    # Update 'rtdataReactive' ..................................................
    rtdataReactive(rtdata)

    # Update feature table
    update_features_ui(session)


    # Render the updated interactive graph
    output$manualFilterPlot_combined <- plotly::renderPlotly({

      rtdata <- rtdataReactive()

      if (!is.null(rtdata$filtering$data)) {
        rtdata$filtering$data <- create_other_isobars(rtdata$filtering$data,
                                                      typeanalysis = "filtering")
      }
      if (!is.null(rtdata$prediction$candidates)) {
        rtdata$prediction$candidates <- create_other_isobars(rtdata$prediction$candidates,
                                                             typeanalysis = "prediction")
      }


      rtdataReactive(rtdata)

      req(lipid_classes_combined())
      req(selected_lipid_classes_combined())

      lipid_classes_combined()
      selected_lipid_classes_combined()

      renderFilterPlot_combined()
    })

    # Notification: Saved successfully (lasts 3 seconds) .......................
    removeNotification(notif_id)  # Remove the 'saving' notification
    showNotification("The changes made in the Highlight have been successfully saved.",
                     type = "message", duration = 3)
  })


  #============================================================================#
  # -- Filtered Combined Interactive Graph
  #============================================================================#
  renderFilterPlot_combined <- function() {

    req(lipid_classes_combined())
    lipidclasses <- lipid_classes_combined()

    rtdata <- rtdataReactive()

    req(selected_lipid_classes_combined())
    lipidclasses_selected <- selected_lipid_classes_combined()
    score <- rtdata$metadata$prediction$rtscorethr
    rt_thr <- as.numeric(input$rt_thr_predict)


    # SHAPES (max 7 subclasses) ................................................
    req(length(lipidclasses) <= 7)
    class_shapes <- c(0, 1, 2, 3, 5, 6, 7)
    shape_map <- setNames(class_shapes[seq_along(lipidclasses)], lipidclasses)


    # Filtering data selected ..................................................
    if ("Filter_final" %in% colnames(rtdata$filtering$data)) {
      rtdata$filtering$data$checkfilter <- rtdata$filtering$data$Filter_final
    }

    filter_data <- reactive({
      rtdata$filtering$data[
        rtdata$filtering$data$LipidSubclass %in% lipidclasses_selected &
          !(rtdata$filtering$data$checkfilter %in% c("Ignore duplicated",
                                                     "Not pass QC", "Unsupported class")),
      ]
    })


    # Prepare model data -------------------------------------------------------

    filter_data_plot <- filter_data()

    filter_data_plot$dataset_type <- "model"

    if (!"other_isobars" %in% names(filter_data_plot)) {
      filter_data_plot$other_isobars <- ""
    }

    if (!"comments" %in% names(filter_data_plot)) {
      filter_data_plot$comments <- ""
    }


    # Prepare prediction data --------------------------------------------------

    prediction_candidates_plot <- NULL

    has_prediction <- (!is.null(rtdata$prediction$candidates) &&
                         nrow(rtdata$prediction$candidates) > 0)

    if (has_prediction) {

      prediction_candidates_plot <- rtdata$prediction$candidates

      if ("Filter_final" %in% colnames(prediction_candidates_plot)) {
        prediction_candidates_plot$checkcandidate <- prediction_candidates_plot$Filter_final
      }

      prediction_candidates_plot <- prediction_candidates_plot[
        prediction_candidates_plot$LipidSubclass %in% lipidclasses_selected &
          !is.na(prediction_candidates_plot$RTsec_rtpred) &
          prediction_candidates_plot$checkcandidate %in%
          c("Kept", "Kept manually by the user"),
      ]


      # prediction_candidates_plot <- prediction_candidates_plot[
      #   !is.na(prediction_candidates_plot$RTsec_rtpred) &
      #     abs(
      #       prediction_candidates_plot$RTsec -
      #         prediction_candidates_plot$RTsec_rtpred
      #     ) <= rt_thr,
      #   ,
      #   drop = FALSE
      # ]

      # Only continue if rows remain after filtering
      if (nrow(prediction_candidates_plot) > 0) {

        prediction_candidates_plot$dataset_type <- "prediction"

        if (!"confidenceLevel" %in% names(prediction_candidates_plot)) {
        prediction_candidates_plot$confidenceLevel <- "Prediction"
        }
        if (!"checkfilter" %in% names(prediction_candidates_plot)) {
        prediction_candidates_plot$checkfilter <- "Kept"
        }
        if (!"Filter" %in% names(prediction_candidates_plot)) {
        prediction_candidates_plot$Filter <- "Prediction"
        }
        if (!"IS" %in% names(prediction_candidates_plot)) {
          prediction_candidates_plot$IS <- ""
        }
        if (!"other_isobars" %in% names(prediction_candidates_plot)) {
          prediction_candidates_plot$other_isobars <- ""
        }
        if (!"comments" %in% names(prediction_candidates_plot)) {
          prediction_candidates_plot$comments <- ""
        }

      } else {

        prediction_candidates_plot <- NULL
      }
    }

    # Combine datasets ---------------------------------------------------------

    common_cols <- c(
      "lipidname", "adduct", "IS", "mz", "RTsec", "C", "DB", "LipidSubclass",
      "confidenceLevel", "Filter", "checkfilter",
      "index", "LipidRTid", "RTsec_rtpred",
      "Score_rtpred", "other_isobars",
      "adducts_expanded", "comments", "dataset_type"
    )

    data_plot <- filter_data_plot[, common_cols]

    if (!is.null(prediction_candidates_plot) &&
        nrow(prediction_candidates_plot) > 0) {

      data_plot <- rbind(
        data_plot,
        prediction_candidates_plot[, common_cols]
      )
    }


    # By default, all DBs have alpha value = 1.
    data_plot$alpha <- ifelse(data_plot$DB %in% input$db_values_filter_combined, 1, 0.8)

    req(!is.null(data_plot))

     req(nrow(data_plot) > 0)
    req(!is.null(input$c_range_filter_combined))
    req(!is.null(input$rt_range_filter_combined))
    req(length(input$c_range_filter_combined) == 2)
    req(length(input$rt_range_filter_combined) == 2)


    # Levels ...................................................................
    kept_levels <- c("Kept",
                     "Kept manually by the user",
                     "Recovered with RT modelling")

    removed_levels <- c("Removed with RT modelling",
                        "In-source fragments",
                        "Removed distant data",
                        "Confidence level among isobars",
                        "Removed by RT rules",
                        "Removed manually by the user",
                        "Outlier",
                        "RT coherence between subclasses",
                        "Discarded - In-source fragments",
                        "Discarded - RT coherence between subclasses",
                        "Discarded - Removed by RT rules",
                        "Discarded - Outlier",
                        "Discarded - Distant data")




    # Range of minimum and maximum values of C from the data. ..................
    # Range of C safely
    if (nrow(data_plot) == 0) return(NULL)

    min_c <- min(data_plot$C, na.rm = TRUE) - 1
    max_c <- max(data_plot$C, na.rm = TRUE) + 1

    if (!is.finite(min_c) || !is.finite(max_c)) return(NULL)

    min_c <- ifelse(min_c %% 2 == 0, min_c - 2, min_c - 1)
    max_c <- ifelse(max_c %% 2 == 0, max_c + 2, max_c + 1)


    minC_filter_combined(min_c)
    maxC_filter_combined(max_c)

    minc_plot <- as.numeric(input$c_range_filter_combined[1])
    maxc_plot <- as.numeric(input$c_range_filter_combined[2])

    req(is.finite(minc_plot), is.finite(maxc_plot))
    req(length(minc_plot) == 1, length(maxc_plot) == 1)


    # Range of minimum and maximum RT values from the data .....................
    min_rt <- floor(min(data_plot$RTsec, na.rm = TRUE))
    max_rt <- ceiling(max(data_plot$RTsec, na.rm = TRUE))

    if (!is.finite(min_rt) || !is.finite(max_rt)) return(NULL)


    minRT_filter_combined(min_rt)
    maxRT_filter_combined(max_rt)

    minrt_plot <- as.numeric(input$rt_range_filter_combined[1])
    maxrt_plot <- as.numeric(input$rt_range_filter_combined[2])

    req(is.finite(minrt_plot), is.finite(maxrt_plot))
    req(length(minrt_plot) == 1, length(maxrt_plot) == 1)



    # Apply graph style --------------------------------------------------------

    if (input$graph_style_filter_combined == "Kept Data") {

      data_plot <- data_plot[
        data_plot$dataset_type == "model" &
          data_plot$checkfilter %in% kept_levels,
      ]

    } else if (input$graph_style_filter_combined == "Removed Data") {

      data_plot <- data_plot[
        data_plot$dataset_type == "model" &
          data_plot$checkfilter %in% removed_levels,
      ]

    } else if (input$graph_style_filter_combined == "Kept Data + Removed Data") {

      data_plot <- data_plot[
        data_plot$dataset_type == "model" &
          data_plot$checkfilter %in%
          c(kept_levels, removed_levels),
      ]

    } else if (input$graph_style_filter_combined == "Kept Data + Accepted Candidates") {

      data_plot <- data_plot[
        (
          data_plot$dataset_type == "model" &
            data_plot$checkfilter %in% kept_levels
        ) |
          (
            data_plot$dataset_type == "prediction"
          ),
      ]

    } else if (input$graph_style_filter_combined == "Kept Data + Removed Data + Accepted Candidates") {

      data_plot <- data_plot[
        (
          data_plot$dataset_type == "model" &
            data_plot$checkfilter %in% c(kept_levels, removed_levels)
        ) |
          (
            data_plot$dataset_type == "prediction"
          ),
      ]

    } else if (input$graph_style_filter_combined == "Accepted Candidates") {

      data_plot <- data_plot[
        data_plot$dataset_type == "prediction",
      ]
    }

    req(nrow(data_plot) > 0)

    if (is.null(data_plot$comments)) {
      data_plot$comments <- ""
    }


    # otrer isobars ............................................................
    # Split por LipidRTid
    # data_plot$other_isobars_text <- vapply(
    #   data_plot$other_isobars,
    #   function(x) {
    #
    #     if (is.na(x) || x == "") {
    #       return("")
    #     }
    #
    #     parts <- strsplit(x, "\\|")[[1]]
    #
    #     lines <- vapply(parts, function(y) {
    #       p <- strsplit(y, ",")[[1]]
    #       paste0(
    #         p[1], " --- ", p[2],
    #         " (", p[3], ")"
    #       )
    #     }, character(1))
    #
    #     paste0(
    #       "<br><b>- Other isobars:</b>",
    #       "<br>   ",
    #       paste(lines, collapse = "<br>   ")
    #     )
    #   },
    #   character(1)
    # )


    # otrer isobars
    # Split por LipidRTid
    # data_plot$other_isobars_text <- vapply(
    #   data_plot$other_isobars,
    #   function(x) {
    #
    #     if (is.na(x) || x == "") {
    #       return("")
    #     }
    #
    #     parts <- strsplit(x, "\\|")[[1]]
    #
    #     lines <- vapply(parts, function(y) {
    #       p <- strsplit(y, ",")[[1]]
    #       paste0(
    #         p[1], " --- ", p[2],
    #         " (", p[3], ")"
    #       )
    #     }, character(1))
    #
    #     paste0(
    #       "<br><b>- Other isobars:</b>",
    #       "<br>   ",
    #       paste(lines, collapse = "<br>   ")
    #     )
    #   },
    #   character(1)
    # )

    # data_plot$other_isobars_text <- vapply(
    #   data_plot$other_isobars,
    #   function(x) {
    #
    #     if (is.na(x) || x == "") {
    #       return("")
    #     }
    #
    #     parts <- strsplit(x, "\\|")[[1]]
    #
    #     lines <- vapply(parts, function(y) {
    #       p <- strsplit(y, ",")[[1]]
    #       paste0(
    #         p[1], " --- ", p[2], "<br>",
    #         " &nbsp;&nbsp; - RTsec: ", p[3], "<br>",
    #         " &nbsp;&nbsp; - RTsec_rtpred: ", p[4], "<br>",
    #         " &nbsp;&nbsp; - Score_rtpred: ", p[5], "<br>"
    #       )
    #     }, character(1))
    #
    #     paste0(
    #       "<br><b>- Other isobars:</b>",
    #       "<br>   ",
    #       paste(lines, collapse = "<br>   ")
    #     )
    #   },
    #   character(1)
    # )
    data_plot$other_isobars_text <- vapply(
      data_plot$other_isobars,
      function(x) {

        if (is.na(x) || x == "") {
          return("")
        }

        parts <- strsplit(x, "\\|")[[1]]

        lines <- vapply(seq_along(parts), function(i) {
          p <- strsplit(parts[i], ",")[[1]]
          paste0(
            "<b> Other Isobar ", i, " -------------------- ", "</b><br>",
            " - ",p[1], " --- ", p[2], "<br>",
            " - RTsec_rtpred: ", p[4], "<br>",
            " - Score_rtpred: ", p[5], "<br>",
            "<b> ----------------------------------------- ", "</b><br>"
          )
        }, character(1))

        paste0(
          "<br>",
          paste(lines, collapse = "<br>")
        )
      },
      character(1)
    )



    # Original Label
    data_plot$original_label_text <- mapply(
      function(filter_value, checkfilter_value) {

        if (!checkfilter_value %in% c(
          "Kept manually by the user",
          "Removed manually by the user",
          "Recovered with RT modelling",
          "Removed with RT modelling"
        )) {
          return("")
        }

        if (is.na(filter_value) || filter_value == "") {
          return("")
        }

        paste0(
          "<br><b>- Original label:</b> ",
          filter_value
        )
      },
      data_plot$Filter,
      data_plot$checkfilter,
      USE.NAMES = FALSE
    )

    data_plot$adduct_text <- vapply(
      data_plot$adducts_expanded,
      function(x) {

        if (is.na(x) || x == "") {
          return("")
        }

        adducts <- strsplit(x, "\\|")[[1]]

        if (length(adducts) == 1) {
          return(
            paste0("<b>", adducts[1], "</b>")
          )
        }

        adducts[1] <- paste0("<b>", adducts[1], "</b>")

        paste(adducts, collapse = " | ")
      },
      character(1)
    )

    # Interactive graphics .....................................................
    lipidname <- IS <- C <- mz <- RTsec <- DB <- index <- LipidRTid <- RTsec_rtpred <- confidenceLevel <- LipidSubclass <- c()
    # If all the DBs are selected...
    if (length(setdiff(sort(unique(filter_data()$DB)), sort(input$db_values_filter_combined))) == 0 &
       !show_filterscore()) {

      p <- ggplot(data = data_plot,
                  aes(x = C, y = RTsec, color = factor(DB), shape = factor(LipidSubclass),
                      text = paste0(
                        "<br><b>", lipidname, "  ---  </b>", adduct_text,
                        ifelse(IS == "IS","<br><b> - Internal Standard</b>",""),

                        # Line
                        "<br>",
                        paste0(
                          rep("─", max(1, nchar(
                            ifelse(IS == "IS",
                                   " - Internal Standard",
                                   paste0(ifelse(is.na(lipidname), "", lipidname), "  ", ifelse(is.na(adduct), "", adduct)))
                          ))),
                          collapse = ""
                        ),
                        "<br><b> - LipidRTid:</b> ", LipidRTid,
                        "<br><b> - mz:</b> ", round(mz,3),
                        "<br><b> - RTsec:</b> ", round(RTsec, 3),
                        "<br><b> - RTsec_rtpred:</b> ", round(RTsec_rtpred, 3),
                        "<br><b> - Score_rtpred:</b> ", Score_rtpred,
                        "<br><b> - confidenceLevel:</b> ", confidenceLevel,
                        "<br><b> - Filter:</b> ", checkfilter,

                        # original label
                        ifelse(
                          !is.na(original_label_text) & original_label_text != "",
                          paste0(
                            "<br>",
                            paste0(
                              rep("─", max(1, nchar(
                                ifelse(IS == "IS",
                                       " - Internal Standard",
                                       paste0(ifelse(is.na(lipidname), "", lipidname),
                                              "  ",
                                              ifelse(is.na(adduct), "", adduct)))
                              ))),
                              collapse = ""
                            ),
                            original_label_text
                          ),
                          ""
                        ),

                        # other isobars
                        ifelse(
                          !is.na(other_isobars_text) & other_isobars_text != "",
                          paste0(
                            "<br>",
                            paste0(
                              rep("─", max(1, nchar(
                                ifelse(IS == "IS",
                                       " - Internal Standard",
                                       paste0(ifelse(is.na(lipidname), "", lipidname),
                                              "  ",
                                              ifelse(is.na(adduct), "", adduct)))
                              ))),
                              collapse = ""
                            ),
                            other_isobars_text
                          ),
                          ""
                        ),

                        # highlight
                        ifelse(
                          !is.na(comments) & comments != "",
                          paste0(
                            "<br>",
                            paste0(
                              rep("─", max(1, nchar(
                                ifelse(IS == "IS",
                                       " - Internal Standard",
                                       paste0(ifelse(is.na(lipidname), "", lipidname),
                                              "  ",
                                              ifelse(is.na(adduct), "", adduct)))
                              ))),
                              collapse = ""
                            ),
                            "<br> [ <b>!</b> ] <b>", comments, "</b> [ <b>!</b> ] "
                          ),
                          ""
                        )
                      ),

                      customdata = index)) +
        geom_point(size = 2.2, stroke = 0.6) +
        labs(x = "Carbon Number", y = "Retention time (s)", shape = "DB") +
        theme_minimal() +
        # theme(legend.position = "none") +
        scale_color_manual(name = "Double \nbonds",
                           values = c("#212121", "#EE2C2C", "#FF9500", "#FFE100",
                                      "#90EE90", "#228B22", "#B0E2FF", "#00B2EE",
                                      "#3A5FCD", "#BF3EFF", "#68228B", "#795548",
                                      "#9E9E9E", "#FFB6C1", "#FFD39B", "#FFF68F",
                                      "#B8860B", "#8B864E", "#B87F9E", "#7FFF00",
                                      "#FF69B4", "#DEBF4B", "#51FFFC", "#EAC4FF",
                                      "#D9FF68",
                                      "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
                                      "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
                                      "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173",
                                      "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
                                      "#969696", "#DD1C77", "#6BAED6", "#74C476", "#9E9AC8",
                                      "#E7969C", "#CEDB9C", "#9C9EDE", "#DBDB8D", "#C49C94",
                                      "#F7B6D2", "#C7C7C7", "#AEC7E8", "#FFBB78", "#98DF8A",
                                      "#FF9896", "#C5B0D5", "#C49C94", "#FDD0A2", "#BDBDBD",
                                      "#5254A3", "#6B6ECF", "#9C9EDE", "#637939", "#B5CF6B",
                                      "#CEDB9C", "#8CA252", "#BD9E39", "#E7BA52", "#E7CB94",
                                      "#843C39", "#AD494A", "#D6616B", "#E7969C", "#7B4173",
                                      "#A55194", "#CE6DBD", "#DE9ED6", "#3182BD", "#6BAED6",
                                      "#9ECAE1", "#C6DBEF", "#E6550D", "#FD8D3C", "#FDAE6B",
                                      "#FDD0A2", "#31A354", "#74C476", "#A1D99B", "#C7E9C0",
                                      "#756BB1", "#9E9AC8", "#BCBDDC", "#DADAEB", "#636363",
                                      "#969696", "#BDBDBD", "#D9D9D9", "#393B79", "#5254A3",
                                      "#6B6ECF", "#9C9EDE", "#637939", "#8CA252", "#B5CF6B",
                                      "#CEDB9C", "#8C6D31", "#BD9E39", "#E7BA52", "#E7CB94"),
                           breaks = seq(0:24)-1) +
        scale_shape_manual(values = shape_map) +
        scale_x_continuous(limits = c(minc_plot, maxc_plot), breaks = seq(minc_plot, maxc_plot, by = 2)) +
        scale_y_continuous(limits = c(minrt_plot, maxrt_plot)) +

        geom_point(
          data = subset(data_plot, !is.na(comments) & comments != ""),
          aes(x = C, y = RTsec),
          shape = 22,
          size = 5,
          fill = "black",
          alpha = 0.25,
          stroke = 0,
          color = NA,
          inherit.aes = FALSE
        ) +

        geom_point(
          data = data_plot[data_plot$dataset_type == "model" &
                             data_plot$checkfilter %in% removed_levels, ],
          aes(x = C, y = RTsec),
          shape = 21,
          size = 5,
          fill = "#EE2C2C",
          alpha = 0.25,
          stroke = 0,
          color = NA,
          inherit.aes = FALSE
        ) +

        geom_point(
          data = data_plot[data_plot$dataset_type == "prediction" &
                             data_plot$checkcandidate %in% removed_levels, ],
          aes(x = C, y = RTsec),
          shape = 21,
          size = 5,
          fill = "#EE2C2C",
          alpha = 0.25,
          stroke = 0,
          color = NA,
          inherit.aes = FALSE
        ) +

        # geom_point(
        #   data = data_plot[data_plot$dataset_type == "prediction", ],
        #   aes(x = C, y = RTsec),
        #   shape = 22,
        #   size = 5,
        #   # fill = "#7f1734",
        #   fill = "#FFD700",
        #   alpha = 0.25,
        #   stroke = 0,
        #   color = NA,
        #   inherit.aes = FALSE
        # ) +

        guides(color = "none", alpha = "none", shape = guide_legend(title = "LipidSubclass"))

    } else { # If any of the DBs are not selected, take alpha into account
      p <- ggplot(data = data_plot,
                  aes(x = C, y = RTsec, color = factor(DB), shape = factor(LipidSubclass),
                      text = paste0(
                        "<br><b>", lipidname, "  ---  </b>", adduct_text,
                        ifelse(IS == "IS","<br><b> - Internal Standard</b>",""),

                        # Line
                        "<br>",
                        paste0(
                          rep("─", max(1, nchar(
                            ifelse(IS == "IS",
                                   " - Internal Standard",
                                   paste0(ifelse(is.na(lipidname), "", lipidname), "  ", ifelse(is.na(adduct), "", adduct)))
                          ))),
                          collapse = ""
                        ),
                        "<br><b> - LipidRTid:</b> ", LipidRTid,
                        "<br><b> - mz:</b> ", round(mz,3),
                        "<br><b> - RTsec:</b> ", round(RTsec, 3),
                        "<br><b> - RTsec_rtpred:</b> ", round(RTsec_rtpred, 3),
                        "<br><b> - Score_rtpred:</b> ", Score_rtpred,
                        "<br><b> - confidenceLevel:</b> ", confidenceLevel,
                        "<br><b> - Filter:</b> ", checkfilter,

                        # original label
                        ifelse(
                          !is.na(original_label_text) & original_label_text != "",
                          paste0(
                            "<br>",
                            paste0(
                              rep("─", max(1, nchar(
                                ifelse(IS == "IS",
                                       " - Internal Standard",
                                       paste0(ifelse(is.na(lipidname), "", lipidname),
                                              "  ",
                                              ifelse(is.na(adduct), "", adduct)))
                              ))),
                              collapse = ""
                            ),
                            original_label_text
                          ),
                          ""
                        ),

                        # other isobars
                        ifelse(
                          !is.na(other_isobars_text) & other_isobars_text != "",
                          paste0(
                            "<br>",
                            paste0(
                              rep("─", max(1, nchar(
                                ifelse(IS == "IS",
                                       " - Internal Standard",
                                       paste0(ifelse(is.na(lipidname), "", lipidname),
                                              "  ",
                                              ifelse(is.na(adduct), "", adduct)))
                              ))),
                              collapse = ""
                            ),
                            other_isobars_text
                          ),
                          ""
                        ),

                        # highlight
                        ifelse(
                          !is.na(comments) & comments != "",
                          paste0(
                            "<br>",
                            paste0(
                              rep("─", max(1, nchar(
                                ifelse(IS == "IS",
                                       " - Internal Standard",
                                       paste0(ifelse(is.na(lipidname), "", lipidname),
                                              "  ",
                                              ifelse(is.na(adduct), "", adduct)))
                              ))),
                              collapse = ""
                            ),
                            "<br> [ <b>!</b> ] <b>", comments, "</b> [ <b>!</b> ] "
                          ),
                          ""
                        )
                      ),

                      alpha = alpha, customdata = index)) +
        geom_point(size = 2.2, stroke = 0.6) +
        labs(x = "Carbon Number", y = "Retention time (s)", shape = "DB") +
        theme_minimal() +
        # theme(legend.position = "top") +
        scale_color_manual(name = "Double \nbonds",
                           values = c("#212121", "#EE2C2C", "#FF9500", "#FFE100",
                                      "#90EE90", "#228B22", "#B0E2FF", "#00B2EE",
                                      "#3A5FCD", "#BF3EFF", "#68228B", "#795548",
                                      "#9E9E9E", "#FFB6C1", "#FFD39B", "#FFF68F",
                                      "#B8860B", "#8B864E", "#B87F9E", "#7FFF00",
                                      "#FF69B4", "#DEBF4B", "#51FFFC", "#EAC4FF",
                                      "#D9FF68",
                                      "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
                                      "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
                                      "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173",
                                      "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
                                      "#969696", "#DD1C77", "#6BAED6", "#74C476", "#9E9AC8",
                                      "#E7969C", "#CEDB9C", "#9C9EDE", "#DBDB8D", "#C49C94",
                                      "#F7B6D2", "#C7C7C7", "#AEC7E8", "#FFBB78", "#98DF8A",
                                      "#FF9896", "#C5B0D5", "#C49C94", "#FDD0A2", "#BDBDBD",
                                      "#5254A3", "#6B6ECF", "#9C9EDE", "#637939", "#B5CF6B",
                                      "#CEDB9C", "#8CA252", "#BD9E39", "#E7BA52", "#E7CB94",
                                      "#843C39", "#AD494A", "#D6616B", "#E7969C", "#7B4173",
                                      "#A55194", "#CE6DBD", "#DE9ED6", "#3182BD", "#6BAED6",
                                      "#9ECAE1", "#C6DBEF", "#E6550D", "#FD8D3C", "#FDAE6B",
                                      "#FDD0A2", "#31A354", "#74C476", "#A1D99B", "#C7E9C0",
                                      "#756BB1", "#9E9AC8", "#BCBDDC", "#DADAEB", "#636363",
                                      "#969696", "#BDBDBD", "#D9D9D9", "#393B79", "#5254A3",
                                      "#6B6ECF", "#9C9EDE", "#637939", "#8CA252", "#B5CF6B",
                                      "#CEDB9C", "#8C6D31", "#BD9E39", "#E7BA52", "#E7CB94"),
                           breaks = seq(0:24)-1) +
        scale_shape_manual(values = shape_map) +
        scale_x_continuous(limits = c(minc_plot, maxc_plot), breaks = seq(minc_plot, maxc_plot, by = 2)) +
        scale_y_continuous(limits = c(minrt_plot, maxrt_plot)) +

        geom_point(
          data = subset(data_plot, !is.na(comments) & comments != ""),
          aes(x = C, y = RTsec),
          shape = 22,
          size = 5,
          fill = "black",
          alpha = 0.25,
          color = NA,
          inherit.aes = FALSE
        ) +

        geom_point(
          data = data_plot[data_plot$dataset_type == "model" &
                             data_plot$checkfilter %in% removed_levels, ],
          aes(x = C, y = RTsec),
          shape = 21,
          size = 5,
          # fill = "#7f1734",
          fill = "#EE2C2C",
          alpha = 0.25,
          color = NA,
          inherit.aes = FALSE
        ) +

        geom_point(
          data = data_plot[data_plot$dataset_type == "prediction" &
                             data_plot$checkcandidate %in% removed_levels, ],
          aes(x = C, y = RTsec),
          shape = 21,
          size = 5,
          fill = "#EE2C2C",
          alpha = 0.25,
          stroke = 0,
          color = NA,
          inherit.aes = FALSE
        ) +

        # geom_point(
        #   data = data_plot[data_plot$dataset_type == "prediction", ],
        #   aes(x = C, y = RTsec),
        #   shape = 21,
        #   size = 5,
        #   # fill = "#7f1734",
        #   fill = "yellow",
        #   alpha = 0.25,
        #   stroke = 0,
        #   color = NA,
        #   inherit.aes = FALSE
        # ) +

          guides(color = "none", alpha = "none", shape = guide_legend(title = "LipidSubclass"))
    }


    plotly::layout(
      plotly::ggplotly(
        p,
        tooltip = c("text"),
        source = "manualFilterPlot_combined_source"
      ),
      legend = list(
        orientation = "h",
        x = 0.5,
        xanchor = "center",
        y = 1.1,

        itemclick = FALSE,
        itemdoubleclick = FALSE
      ),
      dragmode = "lasso"
    )

  }


  #============================================================================#
  # -- Render the initial Manual Filter plot Combined
  #============================================================================#
  output$manualFilterPlot_combined <- plotly::renderPlotly({

    rtdata <- rtdataReactive()

    if (!is.null(rtdata$filtering$data)) {
      rtdata$filtering$data <- create_other_isobars(rtdata$filtering$data,
                                                    typeanalysis = "filtering")
    }
    if (!is.null(rtdata$prediction$candidates)) {
      rtdata$prediction$candidates <- create_other_isobars(rtdata$prediction$candidates,
                                                           typeanalysis = "prediction")
    }

    rtdataReactive(rtdata)

    req(lipid_classes_combined())
    req(selected_lipid_classes_combined())

    lipid_classes_combined()
    selected_lipid_classes_combined()

    renderFilterPlot_combined()
  })






  #============================================================================#
  # -- "Prediction Results" Graphical interface
  #============================================================================#
  update_prediction_ui <- function(session) {
    rtdata <- rtdataReactive()

    # Extract DB values for each Lipid Class.
    db_choices <- DBextraction(rtdata)
    dbChoicesReactive(db_choices)
    dbSelectionsReactive(db_choices)

    # Render the content of the Prediction Results tab .........................
    output$prediction_ui <- renderUI({

      tagList(
        div(
          style = "margin-top: 0px; border-bottom: 1px solid #ddd; padding-bottom: 0px;",
          div(
            class = "no-width-change-tabs",
            fluidRow(
              column(
                width = 6,
                tagList(
                  tags$div(
                    style = "display: flex; justify-content: flex-end; gap: 10px; margin-top: 0px; margin-bottom: 0px;",

                    tags$div(
                      style = "display: flex; align-items: center; justify-content: center; height: 100%;",

                      actionButton(
                        inputId = "manual_prediction_info_btn",
                        label = tags$img(
                          src = "iconohorizontal_buttoninfo.png",
                          style = "height: auto; width: auto; margin-right: 10px; height: 50px;"
                        ),
                        style = "background-color: transparent; border: none;
                        padding: 0; margin-right: 15px; cursor: pointer;"
                      ),

                      # Initially hidden floating menu."
                      tags$div(
                        id = "manual_prediction_dropdown_menu",
                        style = " display: none; position: fixed; top: 50%;
                        left: 50%; transform: translate(-50%, -50%); height: 90vh;
                        width: 80vw; overflow-y: auto; background-color: rgba(255, 255, 255, 1);
                        padding: 20px; font-size: 12px; text-align: right; border: none;
                        border-radius: 5px;  box-shadow: 0 4px 20px rgba(0,0,0,0.3);
                        z-index: 1000;",

                        # Zoom controls.
                        tags$div(
                          style = "position: sticky; top: 0; background-color: rgba(255,255,255,0.95); z-index: 10; padding-bottom: 5px; text-align: right;",
                          actionButton(
                            "zoom_out_dropdown_prediction", "-",
                            style = "width: 25px; height: 25px; padding: 0; font-size: 16px; vertical-align: middle;"
                          ),
                          span(
                            id = "zoom_level_dropdown_prediction", "100%",
                            style = "font-weight: bold; margin: 0 8px; vertical-align: middle;"
                          ),
                          actionButton(
                            "zoom_in_dropdown_prediction", "+",
                            style = "width: 25px; height: 25px; padding: 0; font-size: 16px; vertical-align: middle;"
                          )
                        ),

                        # Manual Prediction Image
                        tags$img(
                          id = "manualPrediction_dropdown",
                          src = "info_ManualPrediction.png",
                          style = "width: 90%; height: auto; display: block; margin: 0 auto;"
                        ),

                        # Script JavaScript for zoom
                        tags$script(HTML("
      var zoomDropdown_prediction = 100;
      var imgDropdown_prediction = document.getElementById('manualPrediction_dropdown');
      var zoomLevelSpanDropdown_prediction = document.getElementById('zoom_level_dropdown_prediction');

      document.getElementById('zoom_in_dropdown_prediction').onclick = function() {
        if (zoomDropdown_prediction < 200) {
          zoomDropdown_prediction += 10;
          imgDropdown_prediction.style.width = zoomDropdown_prediction + '%';
          zoomLevelSpanDropdown_prediction.textContent = zoomDropdown_prediction + '%';
        }
      };
      document.getElementById('zoom_out_dropdown_prediction').onclick = function() {
        if (zoomDropdown_prediction > 30) {
          zoomDropdown_prediction -= 10;
          imgDropdown_prediction.style.width = zoomDropdown_prediction + '%';
          zoomLevelSpanDropdown_prediction.textContent = zoomDropdown_prediction + '%';
        }
      };
    "))
                      ),

                      # Script to show/hide the dropdown when clicking the circle.
                      tags$script(HTML("
    var manualPredictionBtn = document.getElementById('manual_prediction_info_btn');
    var manualPredictionMenu = document.getElementById('manual_prediction_dropdown_menu');

    manualPredictionBtn.onclick = function(event) {
      event.stopPropagation();
      if (manualPredictionMenu.style.display === 'none') {
        manualPredictionMenu.style.display = 'block';
      } else {
        manualPredictionMenu.style.display = 'none';
      }
    };

    document.addEventListener('click', function(event) {
      if (!manualPredictionMenu.contains(event.target) && !manualPredictionBtn.contains(event.target)) {
        manualPredictionMenu.style.display = 'none';
      }
    });
  "))
                    ),
                  )
                )
              )

            ),

            # Tabs for results and manual filtering
            tabsetPanel(
              type = "tabs",
              # Tab: Manual Filter
              tabPanel(
                "Manual Prediction",
                div(
                  style = "margin-top: 40px; margin-bottom: 40px;",
                  fluidRow(
                    column(2,
                           tags$details(
                             style = "margin-bottom: 15px; padding: 8px 12px; border: 1px solid #dcdcdc;
                             border-radius: 6px; background-color: #fafafa;",

                             tags$summary(
                               style = "cursor: pointer; outline: none;",

                               # Title
                               tags$div(
                                 style = "font-size: 14px; font-weight: bold; color: #333; margin-bottom: 4px;",
                                 icon("chart-area"),
                                 HTML("&nbsp;Graph style")
                               ),

                               # Explanation (always visible)
                               tags$div(
                                 style = "font-size: 11px; color: #555; margin-top: 2px; text-align: center;",
                                 "Choose which data layers are displayed."
                               ),

                               # Click hint (below explanation)
                               tags$div(
                                 style = "text-align: center; font-size: 10px; color: #777; margin-top: 5px;",
                                 "click to expand / collapse"
                               )
                             ),

                             ## CONTENT — solo el select
                             div(
                               style = "margin-top: 10px;",

                               selectInput(
                                 inputId = "graph_style",
                                 label   = NULL,
                                 choices = list(
                                   "Combined views" = c(
                                     "Data + Candidates",
                                     "Data + Accepted Candidates"
                                   ),
                                   "Single layers" = c(
                                     "Data",
                                     "Candidates",
                                     "Accepted Candidates",
                                     "Discarded Candidates"
                                   )
                                 ),
                                 selected = "Data + Candidates",
                                 width = "100%"
                               )
                             )
                           ),

                           # RT and C filtering (stacked, collapsible)
                           tags$details(
                             style = "margin-bottom: 15px; border: 1px solid #dcdcdc; border-radius: 6px;
                             background-color: #fafafa; padding: 8px 12px;",

                             ## Header
                             tags$summary(
                               style = "display: block; cursor: pointer; outline: none;",

                               tags$div(
                                 style = "display: flex; align-items: center; font-size: 14px;
                                 font-weight: bold; color: #333; margin-bottom: 2px;",
                                 icon("filter"),
                                 HTML("&nbsp;&nbsp;RT - C Ranges")
                               ),

                               tags$div(
                                 style = "font-size: 11px; color: #555; line-height: 1.3; margin-bottom: 2px; text-align: center;",
                                 HTML(
                                   "Use these controls to modify the Retention Time and Carbon Number
       ranges displayed."
                                 )
                               ),

                               tags$div(
                                 style = "text-align: center; font-size: 10px; color: #777; margin-top: 5px;",
                                 "click to expand / collapse"
                               )
                             ),

                             div(
                               style = "margin-top: 10px;",

                               # RT range
                               div(
                                 style = "padding: 4px; border: 1px dashed #d0d0d0; border-radius: 4px;
                                 background-color: #ffffff; margin-bottom: 10px;",
                                 uiOutput("rt_range_predict_ui")
                               ),

                               # C range
                               div(
                                 style = "padding: 4px; border: 1px dashed #d0d0d0; border-radius: 4px;
                                 background-color: #ffffff;",
                                 uiOutput("c_range_predict_ui")
                               )
                             )
                           ),

                           # DB selector
                           tags$details(
                             id = "db_details_prediction",
                             style = "margin-bottom: 20px; border: 1px solid #dcdcdc;
                             border-radius: 6px; background-color: #fafafa; padding: 8px 12px;",

                             tags$summary(
                               style = "cursor: pointer; outline: none;",

                               tags$div(
                                 style = "font-size: 14px; font-weight: bold; color: #333; margin-bottom: 4px;",
                                 icon("database"),
                                 HTML("&nbsp;DB values")
                               ),

                               tags$div(
                                 style = "font-size: 11px; color: #555; margin-top: 2px; text-align: center;",
                                 HTML("Select the DB values affected by <b>Keep</b> and <b>Remove</b> actions.")
                               ),

                               tags$div(
                                 style = "text-align: center; font-size: 10px; color: #777; margin-top: 5px;",
                                 "click to expand / collapse"
                               )
                             ),

                             div(
                               style = "margin-top: 15px;",
                               uiOutput("db_prediction_content")
                             )
                           ),


                           div(
                             style = "margin-top: 15px; padding: 10px; border: 1px solid #dcdcdc;
                             border-radius: 6px; background-color: #fafafa;",

                             # Title
                             tags$div(
                               style = "font-size: 14px; font-weight: bold; margin-bottom: 4px;",
                               icon("sliders-h"),
                               HTML("&nbsp;Modify predictions")
                             ),

                             # Keep button
                             actionButton(
                               "keepPredictionButton",
                               label = HTML('<span style="font-size:16px; margin-right:6px;">&#x25B2;</span><b style="font-size:14px;">Keep</b>'),
                               style = "width: 100%; margin-top: 6px; margin-bottom: 6px; background-color: #ecf0f1;
                               color: black; border: 1px solid #d0d0d0; padding: 4px 8px;"),

                             # Remove button
                             actionButton(
                               "removePredictionButton",
                               label = HTML('<span style=\"font-size:16px; margin-right:6px;\">&#x25B3;</span><b style=\"font-size:14px;\">Remove</b>'),
                               style = "width: 100%; margin-bottom: 8px; background-color: #ecf0f1;
                               color: black; border: 1px solid #d0d0d0; padding: 4px 8px;"
                             ),

                             # Save button
                             # actionButton(
                             #   "savePredictionButton",
                             #   label = HTML(paste(icon("save"), '<b style="font-size:14px; margin-left:6px;">Save</b>')),
                             #   style = "width: 100%; margin-bottom: 4px; background-color: #ecf0f1;
                             #   color: black; border: 1px solid #d0d0d0; padding: 4px 8px;"
                             # ),
                             uiOutput("savePredictionButtonUI"),

                             # Save hint text
                             tags$div(
                               style = "font-size: 11px; color: #666; margin-top: 2px; text-align: center;",
                               icon("info-circle"),
                               HTML("&nbsp;Changes are only applied after saving")
                             )
                           ),
  #                          div(
  #                            style = "
  #   margin-top: 15px;
  #   padding: 12px;
  #   border: 1px solid #bcdffb;
  #   border-radius: 8px;
  #   background-color: #eef7ff;
  # ",
  #
  #                            # Title
  #                            tags$div(
  #                              style = "font-size: 14px; font-weight: bold; margin-bottom: 6px;",
  #                              icon("check-double"),
  #                              HTML("&nbsp;Final prediction check")
  #                            ),
  #
  #                            # Description
  #                            tags$div(
  #                              style = "font-size: 12px; color: #444; margin-bottom: 10px;",
  #                              HTML(
  #                                paste0(
  #                                  "Once you have finished reviewing the ",
  #                                  "<b>filtering</b> and <b>prediction</b> results, ",
  #                                  "you can run a final joint filtering step to remove ",
  #                                  "remaining inconsistencies and generate the final set ",
  #                                  "of annotations."
  #                                )
  #                              )
  #                            ),
  #
  #                            # Button
  #                            actionButton(
  #                              "checkPredictionsButton",
  #                              label = HTML(
  #                                paste(
  #                                  icon("magic"),
  #                                  '<b style="font-size:14px; margin-left:6px;">
  #       Check predictions
  #       </b>'
  #                                )
  #                              ),
  #                              style = "
  #     width: 100%;
  #     margin-bottom: 8px;
  #     background-color: #d6ebff;
  #     color: #003b6f;
  #     border: 1px solid #9bc7ef;
  #     padding: 6px 8px;
  #     font-weight: bold;
  #   "
  #                            ),
  #
  #                            # Warning box
  #                            tags$div(
  #                              style = "
  #     font-size: 11px;
  #     color: #7a5200;
  #     background-color: #fff4d6;
  #     border: 1px solid #f2d27a;
  #     border-radius: 5px;
  #     padding: 6px;
  #     margin-top: 4px;
  #   ",
  #                              icon("exclamation-triangle"),
  #                              HTML(
  #                                paste0(
  #                                  "&nbsp;<b>Important:</b> After running this step, ",
  #                                  "you may still manually modify filtering or prediction results ",
  #                                  "and save them. However, prediction models will <b>not</b> be ",
  #                                  "retrained and RT predictions will <b>not</b> be recalculated."
  #                                )
  #                              )
  #                            )
  #                          )

                      uiOutput("finalPredictionCheckUI"),

                    ),
                    column(8,

                      ## Interactive Prediction Graph
                      div(
                        class = "full-height",
                        plotly::plotlyOutput("manualPredictionPlot", height = "700px")
                      ),

                      ## Info box under the graph
                      # div(
                      #   style = "margin-top: 12px; padding: 10px 12px; border: 1px solid #dcdcdc;
                      #   border-left: 4px solid #6c757d; border-radius: 4px; background-color: #f8f9fa;
                      #   font-size: 12px; color: #333;",
                      #   icon("info-circle"),
                      #   tags$strong(" Filtering impact: "),
                      #   "If you modify the filtering, the predictions will be recalculated accordingly."
                      # )
                      uiOutput("filterImpactMessage_prediction")
                    ),

                    column(2,
                           # Legend image of the prediction graph. If you click on the image,
                           # it doubles in size. To return it to its original size, tap anywhere on the screen.
                           div(
                             class = "full-width",
                             tags$img(
                               id = "zoomable-image_prediction",
                               src = "legend_manualPrediction.png",
                               width = "90%",
                               style = "cursor: pointer; transition: transform 0.3s ease-in-out; position: relative; top: 0; left: 0;"
                             ),
                             tags$script(HTML("
                                $(document).ready(function() {
                                  $('#zoomable-image_prediction').click(function(event) {
                                    event.stopPropagation();
                                    $(this).toggleClass('zoomed');
                                  });

                                  $(document).click(function() {
                                    $('#zoomable-image_prediction').removeClass('zoomed');
                                  });
                                });
                              ")),
                                                       tags$style(HTML("
                                #zoomable-image_prediction {
                                  max-width: 100%;
                                }

                                #zoomable-image_prediction.zoomed {
                                  transform: scale(2); /* Zoom x2 */
                                  transform-origin: top left;
                                  position: absolute;
                                  top: 0;
                                  left: 0;
                                }
                             "))
                           ),

                           div(
                             style = "margin-top: 20px; padding: 10px; border: 1px solid #dcdcdc;
                             border-radius: 6px; background-color: #fafafa;",

                             # Title
                             tags$div(
                               style = "font-size: 14px; font-weight: bold; margin-bottom: 4px;",
                               icon("highlighter"),
                               HTML("&nbsp;Highlight data")
                             ),

                             # Subtitle / explanation
                             tags$div(
                               style = "font-size: 11px; color: #555; margin-bottom: 8px; text-align: center;",
                               "Mark or Unmark lipids in the plot to highlight them."
                             ),

                             # Mark button
                             actionButton(
                               "markPredictionButton",
                               label = HTML(paste(icon("tag"), '<b style="font-size: 13px; margin-left: 6px;">Mark</b>')),
                               style = "width: 100%; margin-bottom: 6px; background-color: #ecf0f1;
                               color: black; border: 1px solid #d0d0d0; padding: 4px 8px;"),

                             # Unmark button
                             actionButton(
                               "unmarkPredictionButton",
                               label = HTML(paste(icon("times"),'<b style="font-size: 13px; margin-left: 6px;">Unmark</b>')),
                               style = "width: 100%; background-color: #ecf0f1; color: black;
                               border: 1px solid #d0d0d0; padding: 4px 8px; margin-bottom: 4px;"),

                             # Save mark button
                             actionButton(
                               "saveMarkPredictionButton",
                               label = HTML(paste(icon("save"),'<b style="font-size: 13px; margin-left: 6px;">Save</b>')),
                               style = "width: 100%; background-color: #ecf0f1; color: black;
                               border: 1px solid #d0d0d0; padding: 4px 8px; margin-bottom: 4px;"),

                             # Save hint text
                             tags$div(
                               style = "font-size: 11px; color: #666; margin-top: 2px; text-align: center;",
                               icon("info-circle"),
                               HTML("&nbsp;Changes are only applied after saving")
                             )
                           )

                    )
                  )
                )
              ),

              # Tab: Prediction Results
              tabPanel(
                "Prediction Results",
                div(
                  style = "margin-top: 20px;",
                  plotOutput("predictionPlot", height = "600px") # Automatic prediction graphs
                )
              )
            )
          )
        )
      )
    })

    output$filterImpactMessage_prediction <- renderUI({

      rtdata <- rtdataReactive()

      show_message <-
        "Filter_final" %in% colnames(rtdata$prediction$candidates)

      if (show_message) {
        return(NULL)
      }

      div(
        style = "
      margin-top: 12px;
      padding: 10px 12px;
      border: 1px solid #dcdcdc;
      border-left: 4px solid #6c757d;
      border-radius: 4px;
      background-color: #f8f9fa;
      font-size: 12px;
      color: #333;
    ",
        icon("info-circle"),
        tags$strong(" Filtering impact: "),
        "If you modify the filtering, the final prediction check will be reset and predictions will need to be checked again."
      )
    })

    output$savePredictionButtonUI <- renderUI({

      pending_changes <- !is.null(selectedindices_prediction()) &&
        length(selectedindices_prediction()) > 0

      if (pending_changes) {

        actionButton(
          "savePredictionButton",
          label = HTML(
            paste(
              icon("save"),
              '<b style="font-size:14px; margin-left:6px;">Save</b>'
            )
          ),
          style = "
        width: 100%;
        margin-bottom: 4px;
        background-color: #d6ebff;
        color: #003b6f;
        border: 1px solid #9bc7ef;
        padding: 4px 8px;
        font-weight: bold;
      "
        )

      } else {

        actionButton(
          "savePredictionButton",
          label = HTML(
            paste(
              icon("save"),
              '<b style="font-size:14px; margin-left:6px;">Save</b>'
            )
          ),
          style = "
        width: 100%;
        margin-bottom: 4px;
        background-color: #ecf0f1;
        color: black;
        border: 1px solid #d0d0d0;
        padding: 4px 8px;
      "
        )

      }
    })



    output$finalPredictionCheckUI <- renderUI({

      # if ("Filter_final" %in% colnames(rtdataReactive()$prediction$candidates)) {
      #
      #   div(
      #     style = "margin-top: 15px; padding: 12px;  border: 1px solid #c8e6c9;
      #     border-radius: 8px; background-color: #eef9ef;",
      #
      #     tags$div(
      #       style = "font-size: 14px; font-weight: bold; margin-bottom: 6px;",
      #       icon("check-circle"),
      #       HTML("&nbsp;Final prediction check completed")
      #     ),
      #
      #     tags$div(
      #       style = "font-size: 12px; color: #2e7d32; margin-bottom: 10px;",
      #       HTML(
      #         paste0(
      #           "The final prediction check has already been performed. ",
      #           "The current annotations correspond to the final reviewed dataset."
      #         )
      #       )
      #     ),
      #
      #     actionButton(
      #       "redoPredictionCheckButton",
      #       label = HTML(
      #         paste(icon("undo"), '<b style="font-size:14px; margin-left:6px;"> Redo final check</b>')),
      #       style = "width: 100%; background-color: #ecf0f1; color: black;
      #       border: 1px solid #d0d0d0; padding: 4px 8px; margin-bottom: 4px;"
      #     ),
      #
      #     tags$div(
      #       style = "font-size: 11px; color: #666; margin-top: 6px; text-align: center;",
      #       HTML("This will remove the final prediction check results and allow the step to be run again.")
      #     )
      #   )
      #
      # } else {

        div(
          style = "margin-top: 15px; padding: 12px; border: 1px solid #bcdffb;
          border-radius: 8px; background-color: #eef7ff;",

          tags$div(
            style = "font-size: 14px; font-weight: bold; margin-bottom: 6px;",
            icon("check-double"),
            HTML("&nbsp;Final prediction check")
          ),

          tags$div(
            style = "font-size: 12px; color: #444; margin-bottom: 10px;",
            HTML(
              paste0(
                "Once you have finished reviewing the ",
                "<b>filtering</b> and <b>prediction</b> results, ",
                "you can run a final joint filtering step to remove ",
                "remaining inconsistencies and generate the final set ",
                "of annotations."
              )
            )
          ),

          actionButton(
            "checkPredictionsButton",
            label = HTML(
              paste(
                icon("magic"),
                '<b style="font-size:14px; margin-left:6px;">
              Check predictions
            </b>'
              )
            ),
            style = "
  width: 100%;
  margin-bottom: 8px;
  background-color: #ecf0f1;
  color: #003b6f;
  border: 1px solid #003b6f;
  padding: 4px 8px;
"
          ),

          tags$div(
            style = "
          font-size: 11px;
          color: #7a5200;
          background-color: #fff4d6;
          border: 1px solid #f2d27a;
          border-radius: 5px;
          padding: 6px;
          margin-top: 4px;
        ",
            icon("exclamation-triangle"),
            HTML(
              paste0(
                "&nbsp;<b>Important:</b> After running this step, ",
                "you may still manually modify filtering or prediction results ",
                "and save them."
              )
            )
          )
        )
      # }
    })

  #   actionButton(
  #     "checkPredictionsButton",
  #     label = HTML(
  #       paste(
  #         icon("magic"),
  #         '<b style="font-size:14px; margin-left:6px;">
  #       Check predictions
  #     </b>'
  #       )
  #     ),
  #     style = "
  #   width: 100%;
  #   margin-bottom: 8px;
  #   background-color: #ecf0f1;
  #   color: #003b6f;
  #   border: 1px solid #003b6f;
  #   padding: 4px 8px;
  # "
  #   )

    observeEvent(input$checkPredictionsButton, {

      showModal(
        modalDialog(

          title = tagList(
            icon("exclamation-triangle"),
            "Confirm final prediction check"
          ),

          HTML(
            paste0(
              "<p>",
              "This step performs the <b>final joint filtering</b> using both ",
              "filtering and prediction results.",
              "</p>",

              "<p>",
              "Running this step may modify annotations that were previously ",
              "edited manually by the user.",
              "</p>",

              "<p>",
              "<b>Are you sure you want to continue?</b>",
              "</p>"
            )
          ),

          footer = tagList(
            modalButton("Cancel"),

            actionButton(
              "confirmPredictionCheckButton",
              "Run final check",
              class = "btn-primary"
            )
          ),

          easyClose = TRUE
        )
      )

    })


    observeEvent(input$confirmPredictionCheckButton, {

      removeModal()

      # -------------------------
      # Código actual de
      # checkPredictionsButton
      # -------------------------
      notif_id <- showNotification(HTML("Checking predictions and applying final filters... <span class='loader'></span>"),
                                   type = "default",duration = NULL)

      rtdata <- rtdataReactive()

      rtdata <- LipidRT::filterpredictionsdata(rtdata,
                                               distC_predict = input$distC_predict,
                                               minSpecies_predict = input$minSpecies_predict)

      # Changes in filtering .....................................................
      rtdata <- LipidRT::organizertdata(rtdata = rtdata,
                                        typeanalysis = "filtering")

      # Redo the graphs with the changes
      rtdata <- LipidRT::filteringgraph(rtdata = rtdata)

      # Render automatic graphs
      output$filterPlot <- renderPlot({
        req(input$lipid_class_selector)
        lipid_class <- input$lipid_class_selector

        rtdata$filtering$plots[[lipid_class]]
      })

      if (!is.null(rtdata$prediction$candidates)) {
        # Modify the feature matrix with the changes
        rtdata <- LipidRT::organizertdata(rtdata = rtdata,
                                          typeanalysis = "prediction")

        # Redo the graphs with the changes
        rtdata <- LipidRT::predictgraph(rtdata = rtdata,
                                        eval_predict_plot = TRUE)

        # Render automatic graphs
        output$predictionPlot <- renderPlot({
          req(input$lipid_class_selector)
          lipid_class <- input$lipid_class_selector

          rtdata$prediction$plots[[lipid_class]]
        })
      }

      # Update 'rtdataReactive' ..................................................
      rtdataReactive(rtdata)

      # Update feature table
      update_features_ui(session)
      # update_filter_ui(session) # Update the "Filter Results" graphs.
      # update_filter_combined_ui(session) # Update the "Filter Results" graphs.
      # update_prediction_ui(session) # Update the "Prediction Results" graphs


      # Render the updated interactive graph
      output$manualFilterPlot_combined <- plotly::renderPlotly({

        rtdata <- rtdataReactive()

        if (!is.null(rtdata$filtering$data)) {
          rtdata$filtering$data <- create_other_isobars(rtdata$filtering$data,
                                                        typeanalysis = "filtering")
        }
        if (!is.null(rtdata$prediction$candidates)) {
          rtdata$prediction$candidates <- create_other_isobars(rtdata$prediction$candidates,
                                                               typeanalysis = "prediction")
        }


        rtdataReactive(rtdata)

        req(lipid_classes_combined())
        req(selected_lipid_classes_combined())

        lipid_classes_combined()
        selected_lipid_classes_combined()

        renderFilterPlot_combined()
      })

      rtdataReactive(rtdata)

      removeNotification(notif_id)
      showNotification("Final filtering completed successfully.",
                       type = "message", duration = 3)
    })


    #
    # observeEvent(input$checkPredictionsButton, {
    #
    #   notif_id <- showNotification(HTML("Checking predictions and applying final filters... <span class='loader'></span>"),
    #                                type = "default",duration = NULL)
    #
    #   rtdata <- rtdataReactive()
    #
    #   rtdata <- LipidRT::filterpredictionsdata(rtdata)
    #
    #   # Changes in filtering .....................................................
    #   rtdata <- LipidRT::organizertdata(rtdata = rtdata,
    #                                     typeanalysis = "filtering")
    #
    #   # Redo the graphs with the changes
    #   rtdata <- LipidRT::filteringgraph(rtdata = rtdata)
    #
    #   # Render automatic graphs
    #   output$filterPlot <- renderPlot({
    #     req(input$lipid_class_selector)
    #     lipid_class <- input$lipid_class_selector
    #
    #     rtdata$filtering$plots[[lipid_class]]
    #   })
    #
    #   if (!is.null(rtdata$prediction$candidates)) {
    #     # Modify the feature matrix with the changes
    #     rtdata <- LipidRT::organizertdata(rtdata = rtdata,
    #                                       typeanalysis = "prediction")
    #
    #     # Redo the graphs with the changes
    #     rtdata <- LipidRT::predictgraph(rtdata = rtdata,
    #                                     eval_predict_plot = TRUE)
    #
    #     # Render automatic graphs
    #     output$predictionPlot <- renderPlot({
    #       req(input$lipid_class_selector)
    #       lipid_class <- input$lipid_class_selector
    #
    #       rtdata$prediction$plots[[lipid_class]]
    #     })
    #   }
    #
    #   # Update 'rtdataReactive' ..................................................
    #   rtdataReactive(rtdata)
    #
    #   # Update feature table
    #   update_features_ui(session)
    #   # update_filter_ui(session) # Update the "Filter Results" graphs.
    #   # update_filter_combined_ui(session) # Update the "Filter Results" graphs.
    #   # update_prediction_ui(session) # Update the "Prediction Results" graphs
    #
    #
    #   # Render the updated interactive graph
    #   output$manualFilterPlot_combined <- plotly::renderPlotly({
    #
    #     rtdata <- rtdataReactive()
    #
    #     if (!is.null(rtdata$filtering$data)) {
    #       rtdata$filtering$data <- create_other_isobars(rtdata$filtering$data,
    #                                                     typeanalysis = "filtering")
    #     }
    #     if (!is.null(rtdata$prediction$candidates)) {
    #       rtdata$prediction$candidates <- create_other_isobars(rtdata$prediction$candidates,
    #                                                            typeanalysis = "prediction")
    #     }
    #
    #
    #     rtdataReactive(rtdata)
    #
    #     req(lipid_classes_combined())
    #     req(selected_lipid_classes_combined())
    #
    #     lipid_classes_combined()
    #     selected_lipid_classes_combined()
    #
    #     renderFilterPlot_combined()
    #   })
    #
    #   rtdataReactive(rtdata)
    #
    #   removeNotification(notif_id)
    #   showNotification("Final filtering completed successfully.",
    #                    type = "message", duration = 3)
    # })






    output$db_prediction_content <- renderUI({

      req(dbChoicesReactive(), dbSelectionsReactive(), input$lipid_class_selector)

      current_class <- input$lipid_class_selector

      choices_vec <- as.numeric(dbChoicesReactive()$prediction[[current_class]])
      choices_vec <- sort(unique(choices_vec[!is.na(choices_vec)]))

      selected_db_values <- intersect(
        dbSelectionsReactive()$prediction[[current_class]],
        choices_vec
      )

      div(
        checkboxGroupInput(
          "db_values_prediction",
          label    = NULL,
          choices  = choices_vec,
          selected = selected_db_values,
          inline   = FALSE
        ),

        if (length(choices_vec) == 0) {
          tags$div(
            style = "font-size: 11px; color: #777; margin-top: 6px;",
            "No DB values available"
          )
        }
      )
    })

    # Redo Filter predictions
    observeEvent(input$redoPredictionCheckButton, {

      rtdata <- rtdataReactive()

      if ("Filter_final" %in% colnames( rtdata$filtering$data)) {
      rtdata$filtering$data$Filter_final <- NULL
      }
      if ("Filter_final" %in% colnames(rtdata$prediction$candidates)) {
        rtdata$prediction$candidates$Filter_final <- NULL
      }

      rtdataReactive(rtdata)

      showNotification(
        "Final prediction check has been reset. You can run it again at any time.",
        type = "warning",
        duration = 4
      )
    })


    # Observe user changes and update db_selections dynamically ................
    observeEvent(input$db_values_prediction, {

      req(dbSelectionsReactive())
      db_selections <- dbSelectionsReactive()

      req(input$lipid_class_selector)
      current_class <- input$lipid_class_selector

      # new user selection
      db_selections$prediction[[current_class]] <- input$db_values_prediction
      # Save the user's new selection.
      dbSelectionsReactive(db_selections)
    })


    output$predictionPlot <- renderPlot({
      req(input$lipid_class_selector)

      lipid_class <- input$lipid_class_selector

      # Show a message if no class is selected
      validate(
        need(!is.null(lipid_class) && !is.na(lipid_class) && lipid_class != "",
             "No prediction data is available for the selected lipid class. Therefore, RT prediction plots cannot be displayed at this time. Please select a different class or verify your input data.")
      )

      # If there is a class, display the plot
      plot_obj <- rtdata$prediction$plots[[lipid_class]]

      # If the plot is NULL, show a message
      if (is.null(plot_obj)) {
        validate(
          need(FALSE, "No prediction data is available for the selected lipid class. Therefore, RT prediction plots cannot be displayed at this time. Please select a different class or verify your input data.")
        )
      }

      plot_obj
    })


    # Prediction RT Range Interface Values .....................................
    output$rt_range_predict_ui <- renderUI({
      req(minRT_predict(), maxRT_predict())

      rt_min <- minRT_predict()
      rt_max <- maxRT_predict()

      div(
        style = "margin-bottom: 0px;",

        tags$div(
          style = "font-size: 12px; font-weight: bold; margin-bottom: 2px; text-align: center;",
          "Retention Time Range"
        ),

        sliderInput(
          "rt_range_predict",
          label = NULL,
          min = rt_min,
          max = rt_max,
          value = c(rt_min, rt_max),
          step = 10,
          width = "100%"
        )
      )
    })


    # Prediction C Range Interface Values ......................................
    output$c_range_predict_ui <- renderUI({
      req(minC_predict())
      req(maxC_predict())
      c_min <- minC_predict()
      c_max <- maxC_predict()

      div(
        style = "margin-bottom: 0px;",

      tagList(
        tags$div(
          style = "font-size: 12px; font-weight: bold; margin-bottom: 2px; text-align: center;",
          "Carbon Number Range"
        ),
      sliderInput(
        "c_range_predict",
        label = NULL,
        min = c_min,
        max = c_max,
        value = c(c_min, c_max),
        step = 1,
        width = "100%"
      )))
    })
  }


  #============================================================================#
  # -- 'Manual Prediction' Graph Controls
  #============================================================================#

  # REMOVE PREDICTION ..........................................................
  observeEvent(input$removePredictionButton, {

    # Modified lipid indices in the interactive prediction graph
    selected <- plotly::event_data("plotly_selected", source = "manualPredictionPlot_source")

    rtdata <- rtdataReactive()

    # Keep the indices of selected lipids that match the user-selected DBs
    if (!is.null(selected)) {
      new_indices <- as.numeric(unique(selected$customdata))

      indices <- intersect(
        new_indices,
        rtdata$prediction$candidates$index[
          rtdata$prediction$candidates$LipidSubclass == input$lipid_class_selector &
            rtdata$prediction$candidates$DB %in% input$db_values_prediction
        ]
      )


      if ("Filter_final" %in% colnames(rtdata$prediction$candidates)) {

        # Manually assign 'Removed manually by the user' to the checkcandidate field
        rtdata$prediction$candidates$Filter_final[
          rtdata$prediction$candidates$index %in% indices
        ] <- "Removed manually by the user"

        #   indices <- intersect(
        #     new_indices,
        #     rtdata$filtering$data$index[
        #       rtdata$filtering$data$LipidSubclass == input$lipid_class_selector &
        #       rtdata$filtering$data$DB %in% input$db_values_filter
        #     ]
        #   )
        #
        #   idx <- rtdata$filtering$data$index %in% indices
        #
        #   rtdata$filtering$data$checkfilter[idx] <- ifelse(
        #
        #     ## PRIORITY 1: Removed with RT modelling (keep it as is)
        #     rtdata$filtering$data$checkmodel[idx] == "Removed with RT modelling",
        #
        #     "Removed with RT modelling",
        #
        #     ifelse(
        #
        #       ## PRIORITY 2: manual removal cases
        #       rtdata$filtering$data$Filter[idx] %in% c("Kept") |
        #
        #         (rtdata$filtering$data$checkmodel[idx] %in% c("Recovered with RT modelling") &
        #            rtdata$filtering$data$Filter[idx] %in% c(
        #              "Confidence level among isobars",
        #              "Removed by RT rules",
        #              "In-source fragments",
        #              "Removed distant data",
        #              "Outlier",
        #              "RT coherence between subclasses")),
        #
        #       "Removed manually by the user",
        #
        #       ## DEFAULT
        #       rtdata$filtering$data$Filter[idx]
        #     )
        #   )
        #
        #   # Render the updated interactive graph
        #   output$manualFilterPlot <- plotly::renderPlotly({
        #     req(input$lipid_class_selector)
        #     renderFilterPlot()
        #   })
        #

      } else {
        # Manually assign 'Removed manually by the user' to the checkcandidate field
        rtdata$prediction$candidates$checkcandidate[
          rtdata$prediction$candidates$index %in% indices
        ] <-ifelse(
          rtdata$prediction$candidates$checkscore_rtpred[
            rtdata$prediction$candidates$index %in% indices
          ] %in% c("Kept"),
          "Removed manually by the user",
          rtdata$prediction$candidates$checkscore_rtpred[
            rtdata$prediction$candidates$index %in% indices]
        )
      }
    }

    selected_indices <- selectedindices_prediction()
    selected_indices <- c(selected_indices, indices)
    selectedindices_prediction(selected_indices)

    # Update rtdataReactive
    rtdataReactive(rtdata)

    # Render the updated interactive graph
    output$manualPredictionPlot <- plotly::renderPlotly({
      renderPredictionGraph()
    })
  })


  # KEEP PREDICTION ............................................................
  observeEvent(input$keepPredictionButton, {

    # Modified lipid indices in the interactive prediction graph
    selected <- plotly::event_data("plotly_selected", source = "manualPredictionPlot_source")

    rtdata <- rtdataReactive()

    # Keep the indices of selected lipids that match the user-selected DBs
    if (!is.null(selected)) {
      new_indices <- as.numeric(unique(selected$customdata))

      # prediction
      indices <- intersect(
        new_indices,
        rtdata$prediction$candidates$index[
          rtdata$prediction$candidates$LipidSubclass == input$lipid_class_selector &
            rtdata$prediction$candidates$DB %in% input$db_values_prediction
        ]
      )


      if ("Filter_final" %in% colnames(rtdata$prediction$candidates)) {

        # Manually assign 'Kept manually by the user' to the checkcandidate field
        rtdata$prediction$candidates$Filter_final[
          rtdata$prediction$candidates$index %in% indices
        ] <- "Kept manually by the user"
#
#
#         indices <- intersect(
#           new_indices,
#           rtdata$filtering$data$index[
#             rtdata$filtering$data$LipidSubclass == input$lipid_class_selector &
#             rtdata$filtering$data$DB %in% input$db_values_filter
#           ]
#         )
#
#         idx <- rtdata$filtering$data$index %in% indices
#
#         checkmodel <- rtdata$filtering$data$checkmodel[idx]
#         Filter     <- rtdata$filtering$data$Filter[idx]
#
#         rtdata$filtering$data$checkfilter[idx] <- ifelse(
#
#           ## PRIORITY 1: keep "Recovered with RT modelling"
#           checkmodel == "Recovered with RT modelling",
#           "Recovered with RT modelling",
#
#           ifelse(
#
#             ## PRIORITY 2: manual keep cases
#             Filter %in% c("Confidence level among isobars", "Removed by RT rules",
#                           "In-source fragments", "Removed distant data",
#                           "Outlier", "RT coherence between subclasses") |
#               (checkmodel == "Removed with RT modelling" & Filter == "Kept"),
#
#             "Kept manually by the user",
#
#             ## DEFAULT
#             Filter
#           )
#         )
#
#       # Render the updated interactive graph
#       output$manualFilterPlot <- plotly::renderPlotly({
#         req(input$lipid_class_selector)
#         renderFilterPlot()
#       })

      } else {
        # Manually assign 'Kept manually by the user' to the checkcandidate field
        rtdata$prediction$candidates$checkcandidate[
          rtdata$prediction$candidates$index %in% indices
        ] <- ifelse(
          rtdata$prediction$candidates$checkscore_rtpred[
            rtdata$prediction$candidates$index %in% indices
          ] %in% c("Discarded", "Discarded - In-source fragments",
                   "Discarded - RT coherence between subclasses",
                   "Discarded - Removed by RT rules", "Discarded - Outlier",
                   "Discarded - Distant data"),
          "Kept manually by the user",
          rtdata$prediction$candidates$checkscore_rtpred[
            rtdata$prediction$candidates$index %in% indices]
        )

      }
    }

    selected_indices <- selectedindices_prediction()
    selected_indices <- c(selected_indices, indices)
    selectedindices_prediction(selected_indices)

    # Update 'rtdataReactive'
    rtdataReactive(rtdata)

    # Render the updated interactive graph
    output$manualPredictionPlot <- plotly::renderPlotly({
      renderPredictionGraph()
    })
  })


  # MARK PREDICTION ............................................................
  # observeEvent(input$markPredictionButton, {
  #
  #   # Modified lipid indices in the interactive prediction graph
  #   selected <- plotly::event_data("plotly_selected", source = "manualPredictionPlot_source")
  #
  #   rtdata <- rtdataReactive()
  #
  #   # Keep the indices of selected lipids that match the user-selected DBs
  #   if (!is.null(selected)) {
  #     new_indices <- as.numeric(unique(selected$customdata))
  #
  #     indices <- intersect(
  #       new_indices,
  #       rtdata$prediction$candidates$index[
  #         rtdata$prediction$candidates$LipidSubclass == input$lipid_class_selector &
  #           rtdata$prediction$candidates$DB %in% input$db_values_prediction
  #       ]
  #     )
  #
  #     # Manually assign 'Kept manually by the user' to the checkcandidate field
  #     rtdata$prediction$candidates$comments[
  #       rtdata$prediction$candidates$index %in% indices
  #     ] <- "Marked"
  #   }
  #
  #   # Update 'rtdataReactive'
  #   rtdataReactive(rtdata)
  #
  #   # Render the updated interactive graph
  #   output$manualPredictionPlot <- plotly::renderPlotly({
  #     renderPredictionGraph()
  #   })
  # })
  observeEvent(input$markPredictionButton, {

    showModal(modalDialog(
      title = NULL,

      div(
        style = "
        padding: 15px;
        background-color: #f8f9fa;
        border-radius: 8px;
      ",

        h4(icon("tag"), " Mark selected candidates"),
        p(
          style = "color: #666;",
          "Enter a comment that will be attached to all selected candidates."
        ),

        textAreaInput(
          "mark_text",
          label = NULL,
          placeholder = "Write your comment here...",
          rows = 3
        )
      ),

      footer = tagList(
        modalButton("Cancel"),
        actionButton(
          "confirm_mark_prediction",
          "Apply comment",
          class = "btn-primary"
        )
      ),

      easyClose = TRUE
    ))
  })

  observe({
    shinyjs::toggleState(
      "confirm_mark_prediction",
      condition = nzchar(trimws(input$mark_text))
    )
  })


  observeEvent(input$confirm_mark_prediction, {

    removeModal()

    # Modified lipid indices in the interactive prediction graph
    selected <- plotly::event_data("plotly_selected", source = "manualPredictionPlot_source")

    rtdata <- rtdataReactive()

    # Keep the indices of selected lipids that match the user-selected DBs
    if (!is.null(selected)) {
      new_indices <- as.numeric(unique(selected$customdata))

      indices <- intersect(
        new_indices,
        rtdata$prediction$candidates$index[
          rtdata$prediction$candidates$LipidSubclass == input$lipid_class_selector &
            rtdata$prediction$candidates$DB %in% input$db_values_prediction
        ]
      )

      # Manually assign 'Kept manually by the user' to the checkcandidate field
      rtdata$prediction$candidates$comments[
        rtdata$prediction$candidates$index %in% indices
      ] <- input$mark_text
    }

    # Update 'rtdataReactive'
    rtdataReactive(rtdata)

    # Render the updated interactive graph
    output$manualPredictionPlot <- plotly::renderPlotly({
      renderPredictionGraph()
    })
  })

  # UNMARK PREDICTION ..........................................................
  observeEvent(input$unmarkPredictionButton, {

    # Modified lipid indices in the interactive prediction graph
    selected <- plotly::event_data("plotly_selected", source = "manualPredictionPlot_source")

    rtdata <- rtdataReactive()

    # Keep the indices of selected lipids that match the user-selected DBs
    if (!is.null(selected)) {
      new_indices <- as.numeric(unique(selected$customdata))

      indices <- intersect(
        new_indices,
        rtdata$prediction$candidates$index[
          rtdata$prediction$candidates$LipidSubclass == input$lipid_class_selector &
            rtdata$prediction$candidates$DB %in% input$db_values_prediction
        ]
      )

      rtdata$prediction$candidates$comments[
        rtdata$prediction$candidates$index %in% indices
      ] <- ""
    }

    # Update 'rtdataReactive'
    rtdataReactive(rtdata)

    # Render the updated interactive graph
    output$manualPredictionPlot <- plotly::renderPlotly({
      renderPredictionGraph()
    })
  })


  # SAVE MARK PREDICTION ............................................................
  observeEvent(input$saveMarkPredictionButton, {
    notif_id <- showNotification(HTML("The changes from the Highlight are being saved... <span class='loader'></span>"),
                                 type = "default", duration = NULL)

    # Save the changes to the reactive rtdata.
    rtdata <- rtdataReactive()

    # Changes in filtering .....................................................
    # Modify the feature matrix with the changes
    rtdata <- LipidRT::organizertdata(rtdata = rtdata,
                                      typeanalysis = "prediction")

    # Redo the graphs with the changes
    rtdata <- LipidRT::predictgraph(rtdata = rtdata,
                                    eval_predict_plot = TRUE)

    # Update 'rtdataReactive'
    rtdataReactive(rtdata)

    # Update feature table
    update_features_ui(session)

    # Render automatic graphs
    output$predictionPlot <- renderPlot({
      req(input$lipid_class_selector)
      lipid_class <- input$lipid_class_selector

      rtdata$prediction$plots[[lipid_class]]
    })

    # Render the updated interactive graph
    output$manualPredictionPlot <- plotly::renderPlotly({
      req(input$lipid_class_selector)

      renderPredictionGraph()
    })

    removeNotification(notif_id)  # Remove the 'saving' notification.
    showNotification("The changes made in the Highlight have been successfully saved.", type = "message", duration = 2)

  })




  # SAVE PREDICTION ............................................................
  observeEvent(input$savePredictionButton, {
    notif_id <- showNotification(HTML("The changes from the predictions are being saved... <span class='loader'></span>"),
                                 type = "default", duration = NULL)

    # Save the changes to the reactive rtdata.
    rtdata <- rtdataReactive()

    # Changes in filtering .....................................................

    # Filter prediction results  - no hacer checkpredictions
    # rtdata <- filterpredictionsdata(rtdata,
    #                                 predict_plot = FALSE)
    #
    # rtdata$filtering$data <- create_other_isobars(rtdata$filtering$data,
    #                                               typeanalysis = "filtering")
    #
    # # Render automatic graphs
    # output$filterPlot <- renderPlot({
    #   req(input$lipid_class_selector)
    #   lipid_class <- input$lipid_class_selector
    #
    #   rtdata$filtering$plots[[lipid_class]]
    # })
    #
    # # Render the updated interactive graph
    # output$manualFilterPlot <- plotly::renderPlotly({
    #   req(input$lipid_class_selector)
    #
    #   renderFilterPlot()
    # })

    rtdata <- predictgraph(rtdata = rtdata,
                           classes = NULL,
                           eval_predict_plot = TRUE)

    rtdata <- organizertdata(rtdata = rtdata,
                             typeanalysis = "prediction")

    rtdata$filtering$data <- create_other_isobars(rtdata$filtering$data,
                                                  typeanalysis = "filtering")

    rtdata$prediction$candidates <- create_other_isobars(rtdata$prediction$candidates,
                                                         typeanalysis = "prediction")

    # Update feature table
    update_features_ui(session)

    # Render automatic graphs
    output$predictionPlot <- renderPlot({
      req(input$lipid_class_selector)
      lipid_class <- input$lipid_class_selector

        rtdata$prediction$plots[[lipid_class]]
    })

    # Render the updated interactive graph
    output$manualPredictionPlot <- plotly::renderPlotly({
      req(input$lipid_class_selector)

      renderPredictionGraph()
    })


    # Update 'rtdataReactive'
    rtdataReactive(rtdata)


    output$manualFilterPlot_combined <- plotly::renderPlotly({


      req(lipid_classes_combined())
      req(selected_lipid_classes_combined())

      lipid_classes_combined()
      selected_lipid_classes_combined()

      renderFilterPlot_combined()
    })

    selectedindices_prediction(NULL)

    # Update feature table
    update_features_ui(session)


    removeNotification(notif_id)  # Remove the 'saving' notification.

    if ("Filter_final" %in% colnames(rtdata$prediction$candidates)) {
      showNotification("The changes made have been successfully saved.", type = "message", duration = 2)
    } else {
      showNotification("The changes made in the predictions have been successfully saved.", type = "message", duration = 2)
    }

  })


  #============================================================================#
  # -- Prediction Interactive Graph
  #============================================================================#
  renderPredictionGraph <- function() {

    req(rtdataReactive())
    req(input$lipid_class_selector)

    rtdata <- rtdataReactive()
    lipidclass <- input$lipid_class_selector
    score <- rtdata$metadata$prediction$rtscorethr
    rt_thr <- as.numeric(input$rt_thr_predict)

    # Prediction data ..........................................................
    if ("Filter_final" %in% colnames(rtdata$filtering$data)) {
      rtdata$filtering$data$checkfilter <-  rtdata$filtering$data$Filter_final
    }

    prediction_data <- reactive({
        rtdata$filtering$data[
          rtdata$filtering$data$LipidSubclass == lipidclass,
        ]
    })

    # Get the candidates .......................................................
    if ("Filter_final" %in% colnames(rtdata$prediction$candidates)) {
      rtdata$prediction$candidates$checkcandidate <-  rtdata$prediction$candidates$Filter_final
    }

      candidates_data <- reactive({
        rtdata$prediction$candidates[
          rtdata$prediction$candidates$LipidSubclass == lipidclass &
            !is.na(rtdata$prediction$candidates$RTsec_rtpred) &
            rtdata$prediction$candidates$checkcandidate %in%
            c("Kept","Discarded", "Discarded - In-source fragments",
              "Discarded - RT coherence between subclasses",
              "Discarded - Removed by RT rules", "Discarded - Outlier",
              "Discarded - Distant data",
              "Kept manually by the user","Removed manually by the user"),
        ]
      })


      cols <- c(
        "lipidname", "adduct", "IS", "mz", "RTsec", "C", "DB",
        "confidenceLevel", "checkcandidate", "index", "LipidRTid",
        "RTsec_rtpred", "CI_lower", "CI_upper", "Score_rtpred",
        "shape", "fill", "legend_label", "other_isobars",
        "adducts_expanded", "comments"
      )

    # Get the predictions of the candidates ....................................
    # predictions_data <- reactive({
    #     rtdata$prediction$predictions[
    #       rtdata$prediction$predictions$LipidSubclass == lipidclass,
    #     ]
    # })
    # updated_predictions <- predictions_data()

    # Prepare the data for plotting ............................................
    prediction_data_plot <- prediction_data()[
      prediction_data()$checkfilter %in%
        c("Kept",
          "Kept manually by the user",
          "Recovered with RT modelling"),, drop = FALSE]

    if (!is.null(prediction_data_plot) && nrow(prediction_data_plot) > 0) {

      prediction_data_plot$checkcandidate <- "Data Model"
      prediction_data_plot$shape <- 16
      prediction_data_plot$fill <- "black"
      prediction_data_plot$legend_label <- "Data Model"
      prediction_data_plot$other_isobars <- ""
      prediction_data_plot$comments <- ""
    } else {
      prediction_data_plot <- setNames(
        data.frame(matrix(ncol = length(cols), nrow = 0)),
        cols
      )
    }

    prediction_candidates_plot <- candidates_data()

    # If there are candidates and predictions have been made...
    if (!is.null(prediction_candidates_plot) && nrow(prediction_candidates_plot) > 0) {
      prediction_candidates_plot$shape <- 22
      prediction_candidates_plot$fill <- ifelse(prediction_candidates_plot$checkcandidate %in%
                                                  c("Kept", "Kept manually by the user"),
                                                "white", "grey")
      prediction_candidates_plot$legend_label <- ifelse(prediction_candidates_plot$checkcandidate %in%
                                                          c("Kept", "Kept manually by the user"),
                                                        "Accepted Candidates",
                                                        "Discarded Candidates")


        # prediction_candidates_plot <- prediction_candidates_plot[
        #   !is.na(prediction_candidates_plot$RTsec_rtpred) &
        #     (prediction_candidates_plot$legend_label == "Accepted Candidates" |
        #        (prediction_candidates_plot$legend_label == "Discarded Candidates" &
        #           abs(prediction_candidates_plot$RTsec - prediction_candidates_plot$RTsec_rtpred) <= rt_thr) |
        #        prediction_candidates_plot$checkcandidate %in% c("Kept manually by the user",
        #                                                         "Removed manually by the user")),]
      prediction_candidates_plot <- prediction_candidates_plot[
        !is.na(prediction_candidates_plot$RTsec_rtpred) &
          (prediction_candidates_plot$legend_label %in% c("Accepted Candidates", "Discarded Candidates") |
             prediction_candidates_plot$checkcandidate %in% c("Kept manually by the user",
                                                              "Removed manually by the user")),]



      prediction_candidates_plot$confidenceLevel <- "Prediction"
      prediction_candidates_plot$IS <- ""

      if (!"other_isobars" %in% names(prediction_candidates_plot)) {
        prediction_candidates_plot$other_isobars <- ""
      }
      if (!"comments" %in% names(prediction_candidates_plot)) {
        prediction_candidates_plot$comments <- ""
      }

      prediction_plot <- rbind(
        prediction_data_plot[, cols],
        prediction_candidates_plot[, cols]
      )

    } else {
      # If there are no candidates or no predictions have been made, keep only the model data
      prediction_plot <- prediction_data_plot
    }

    # Convert legend_label to a factor with defined levels.
    prediction_plot$legend_label <- factor(prediction_plot$legend_label,
                                           levels = c("Data Model",
                                                      "Accepted Candidates",
                                                      "Discarded Candidates"))

    # By default, all DBs have alpha value = 1.
    prediction_plot$alpha <- ifelse(prediction_plot$DB %in% input$db_values_prediction, 1, 0.8)

    req(!is.null(prediction_plot))
    req(nrow(prediction_plot) > 0)

    req(!is.null(input$c_range_predict))
    req(!is.null(input$rt_range_predict))
    req(length(input$c_range_predict) == 2)
    req(length(input$rt_range_predict) == 2)

    # C Range ..................................................................
    min_c <- min(prediction_plot$C, na.rm = TRUE) - 1
    max_c <- max(prediction_plot$C, na.rm = TRUE) + 1

    if (!is.finite(min_c) || !is.finite(max_c)) return(NULL)

    min_c <- ifelse(min_c %% 2 == 0, min_c - 2, min_c - 1)
    max_c <- ifelse(max_c %% 2 == 0, max_c + 2, max_c + 1)

    minC_predict(min_c)
    maxC_predict(max_c)

    minc_plot <- as.numeric(input$c_range_predict[1])
    maxc_plot <- as.numeric(input$c_range_predict[2])

    req(is.finite(minc_plot), is.finite(maxc_plot))
    req(length(minc_plot) == 1, length(maxc_plot) == 1)


    # RT Range .................................................................
    min_rt <- floor(min(prediction_plot$RTsec, na.rm = TRUE)) - 50
    max_rt <- ceiling(max(prediction_plot$RTsec, na.rm = TRUE)) + 50

    if (!is.finite(min_rt) || !is.finite(max_rt)) return(NULL)


    minRT_predict(min_rt)
    maxRT_predict(max_rt)

    minrt_plot <- as.numeric(input$rt_range_predict[1])
    maxrt_plot <- as.numeric(input$rt_range_predict[2])

    req(is.finite(minrt_plot), is.finite(maxrt_plot))
    req(length(minrt_plot) == 1, length(maxrt_plot) == 1)

    # prediction lines .........................................................
    # check if prediction layer is required
    prediction_layer <- LipidRT::predictionlayergraph(rtdata = rtdata,
                                                      lipidsubclass = lipidclass,
                                                      min_c = minc_plot,
                                                      max_c = maxc_plot,
                                                      min_rt = minrt_plot,
                                                      max_rt = maxrt_plot,
                                                      alpha = 0.7,
                                                      typeanalysis = "prediction")

    # Define shapes and fills globally .........................................
    shapes <- c("Accepted Candidates" = 22,
                "Discarded Candidates" = 22)
    fills  <- c("Accepted Candidates" = "white",
                "Discarded Candidates" = "grey")

    shapes_present = shapes[names(shapes) %in% unique(prediction_plot$legend_label)]
    fills_present = fills[names(fills) %in% unique(prediction_plot$legend_label)]

    # selected graph style .....................................................
    if (input$graph_style == "Data + Candidates") {
    } else if (input$graph_style == "Data + Accepted Candidates") {
      prediction_plot <- prediction_plot[prediction_plot$legend_label %in% c("Data Model",
                                                                 "Accepted Candidates"),, drop = FALSE]
    } else if (input$graph_style == "Data") {
      prediction_plot <- prediction_plot[prediction_plot$legend_label == "Data Model",, drop = FALSE]
    } else if (input$graph_style == "Candidates") {
      prediction_plot <- prediction_plot[prediction_plot$legend_label %in% c("Accepted Candidates",
                                                                 "Discarded Candidates"),, drop = FALSE]
    } else if (input$graph_style == "Accepted Candidates") {
      prediction_plot <- prediction_plot[prediction_plot$legend_label %in% c("Accepted Candidates"),, drop = FALSE]
    } else if (input$graph_style == "Discarded Candidates") {
      prediction_plot <- prediction_plot[prediction_plot$legend_label %in% c("Discarded Candidates"),, drop = FALSE]
    }

    if (nrow(prediction_plot) == 0) return(NULL)

    # otrer isobars
    # Split por LipidRTid
    # prediction_plot$other_isobars_text <- vapply(
    #   prediction_plot$other_isobars,
    #   function(x) {
    #
    #     if (is.na(x) || x == "") {
    #       return("")
    #     }
    #
    #     parts <- strsplit(x, "\\|")[[1]]
    #
    #     lines <- vapply(parts, function(y) {
    #       p <- strsplit(y, ",")[[1]]
    #       paste0(
    #         p[1], " --- ", p[2],
    #         " (", p[3], ")"
    #       )
    #     }, character(1))
    #
    #     paste0(
    #       "<br><b>- Other isobars:</b>",
    #       "<br>   ",
    #       paste(lines, collapse = "<br>   ")
    #     )
    #   },
    #   character(1)
    # )

    # prediction_plot$other_isobars_text <- vapply(
    #   prediction_plot$other_isobars,
    #   function(x) {
    #
    #     if (is.na(x) || x == "") {
    #       return("")
    #     }
    #
    #     parts <- strsplit(x, "\\|")[[1]]
    #
    #     lines <- vapply(parts, function(y) {
    #       p <- strsplit(y, ",")[[1]]
    #       paste0(
    #         p[1], " --- ", p[2], "<br>",
    #         " &nbsp;&nbsp; - RTsec: ", p[3], "<br>",
    #         " &nbsp;&nbsp; - RTsec_rtpred: ", p[4], "<br>",
    #         " &nbsp;&nbsp; - Score_rtpred: ", p[5], "<br>"
    #       )
    #     }, character(1))
    #
    #     paste0(
    #       "<br><b>- Other isobars:</b>",
    #       "<br>   ",
    #       paste(lines, collapse = "<br>   ")
    #     )
    #   },
    #   character(1)
    # )

    prediction_plot$other_isobars_text <- vapply(
      prediction_plot$other_isobars,
      function(x) {

        if (is.na(x) || x == "") {
          return("")
        }

        parts <- strsplit(x, "\\|")[[1]]

        lines <- vapply(seq_along(parts), function(i) {
          p <- strsplit(parts[i], ",")[[1]]
          paste0(
            "<b> Other Isobar ", i, " -------------------- ", "</b><br>",
            " - ",p[1], " --- ", p[2], "<br>",
            " - RTsec_rtpred: ", p[4], "<br>",
            " - Score_rtpred: ", p[5], "<br>",
            "<b> ----------------------------------------- ", "</b><br>"
          )
        }, character(1))

        paste0(
          "<br>",
          paste(lines, collapse = "<br>")
        )
      },
      character(1)
    )



    prediction_plot$adduct_text <- vapply(
      prediction_plot$adducts_expanded,
      function(x) {

        if (is.na(x) || x == "") {
          return("")
        }

        adducts <- strsplit(x, "\\|")[[1]]

        if (length(adducts) == 1) {
          return(
            paste0("<b>", adducts[1], "</b>")
          )
        }

        adducts[1] <- paste0("<b>", adducts[1], "</b>")

        paste(adducts, collapse = " | ")
      },
      character(1)
    )


    # Initialize variables to avoid errors in aes().
    lipidname <- C <- IS <- mz <- RTsec <- DB <- index <- LipidRTid <- RTsec_rtpred <- CI_lower <- CI_upper <- c()
      Score_rtpred <- confidenceLevel <- checkcandidate <- legend_label <- alpha <- c()
    # Check that prediction_plot has data.
    if (!is.null(prediction_plot) && nrow(prediction_plot) > 0) {
      prediction_plot$C <- as.numeric(as.character(prediction_plot$C))
      # Check that all DBs are in input$db_values_prediction.
      if (length(intersect(unique(prediction_plot$DB), input$db_values_prediction)) == length(unique(prediction_plot$DB))) {

        p <- ggplot(prediction_plot, aes(x = C, y = RTsec, color = as.factor(DB),
                                         text = paste0(
                                           "<br><b>", lipidname, "  ---  </b>", adduct_text,

                                           ifelse(IS == "IS","<br><b> - Internal Standard</b>",""),

                                           # Line
                                           "<br>",
                                           paste0(
                                             rep("─", max(1, nchar(
                                               ifelse(IS == "IS",
                                                      " - Internal Standard",
                                                      paste0(ifelse(is.na(lipidname), "", lipidname), "  ", ifelse(is.na(adduct), "", adduct)))
                                             ))),
                                             collapse = ""
                                           ),
                                           "<br><b> - LipidRTid:</b> ", LipidRTid,
                                           "<br><b> - mz:</b> ", round(mz,3),
                                           "<br><b> - RTsec:</b> ", round(RTsec, 3),
                                           "<br><b> - RTsec_rtpred:</b> ", round(RTsec_rtpred, 3),
                                           "<br><b> - Score_rtpred:</b> ", Score_rtpred,
                                           "<br><b> - confidenceLevel:</b> ", confidenceLevel,
                                           "<br>   (", checkcandidate, ")",

                                           # other isobars
                                           ifelse(
                                             !is.na(other_isobars_text) & other_isobars_text != "",
                                             paste0(
                                               "<br>",
                                               paste0(
                                                 rep("─", max(1, nchar(
                                                   ifelse(IS == "IS",
                                                          " - Internal Standard",
                                                          paste0(ifelse(is.na(lipidname), "", lipidname),
                                                                 "  ",
                                                                 ifelse(is.na(adduct), "", adduct)))
                                                 ))),
                                                 collapse = ""
                                               ),
                                               other_isobars_text
                                             ),
                                             ""
                                           ),

                                           # highlight
                                           ifelse(
                                             !is.na(comments) & comments != "",
                                             paste0(
                                               "<br>",
                                               paste0(
                                                 rep("─", max(1, nchar(
                                                   ifelse(IS == "IS",
                                                          " - Internal Standard",
                                                          paste0(ifelse(is.na(lipidname), "", lipidname),
                                                                 "  ",
                                                                 ifelse(is.na(adduct), "", adduct)))
                                                 ))),
                                                 collapse = ""
                                               ),
                                               "<br> [ <b>!</b> ] <b>", comments, "</b> [ <b>!</b> ] "
                                             ),
                                             ""
                                           )
                                         ),

                                         customdata = index)) +
          geom_point(data = prediction_plot[prediction_plot$shape == 16,], size = 2.2, stroke = 0.6) +
          geom_point(data = prediction_plot[prediction_plot$shape == 22,],
                     aes(shape = legend_label, fill = legend_label), size = 2.2, stroke = 0.6) +
          labs(x = "Carbon number", y = "Retention time (s)") +
          scale_color_manual(name = "Double \nbonds",
                             values = c("#212121", "#EE2C2C", "#FF9500", "#FFE100",
                                        "#90EE90", "#228B22", "#B0E2FF", "#00B2EE",
                                        "#3A5FCD", "#BF3EFF", "#68228B", "#795548",
                                        "#9E9E9E", "#FFB6C1", "#FFD39B", "#FFF68F",
                                        "#B8860B", "#8B864E", "#B87F9E", "#7FFF00",
                                        "#FF69B4", "#DEBF4B", "#51FFFC", "#EAC4FF",
                                        "#D9FF68",
                                        "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
                                        "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
                                        "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173",
                                        "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
                                        "#969696", "#DD1C77", "#6BAED6", "#74C476", "#9E9AC8",
                                        "#E7969C", "#CEDB9C", "#9C9EDE", "#DBDB8D", "#C49C94",
                                        "#F7B6D2", "#C7C7C7", "#AEC7E8", "#FFBB78", "#98DF8A",
                                        "#FF9896", "#C5B0D5", "#C49C94", "#FDD0A2", "#BDBDBD",
                                        "#5254A3", "#6B6ECF", "#9C9EDE", "#637939", "#B5CF6B",
                                        "#CEDB9C", "#8CA252", "#BD9E39", "#E7BA52", "#E7CB94",
                                        "#843C39", "#AD494A", "#D6616B", "#E7969C", "#7B4173",
                                        "#A55194", "#CE6DBD", "#DE9ED6", "#3182BD", "#6BAED6",
                                        "#9ECAE1", "#C6DBEF", "#E6550D", "#FD8D3C", "#FDAE6B",
                                        "#FDD0A2", "#31A354", "#74C476", "#A1D99B", "#C7E9C0",
                                        "#756BB1", "#9E9AC8", "#BCBDDC", "#DADAEB", "#636363",
                                        "#969696", "#BDBDBD", "#D9D9D9", "#393B79", "#5254A3",
                                        "#6B6ECF", "#9C9EDE", "#637939", "#8CA252", "#B5CF6B",
                                        "#CEDB9C", "#8C6D31", "#BD9E39", "#E7BA52", "#E7CB94"),
                             breaks = seq(0:24)-1) +
          scale_x_continuous(limits = c(minc_plot, maxc_plot), breaks = seq(minc_plot, maxc_plot, by = 2)) +
          scale_y_continuous(limits = c(minrt_plot, maxrt_plot)) +
          theme_minimal() +
          theme(legend.position = "none") +

          geom_point(data = subset(prediction_plot, !is.na(comments) & comments != ""),
                     aes(x = C, y = RTsec), shape = 22, size = 5,
                     fill = "black", alpha = 0.25, color = NA, inherit.aes = FALSE)

        # Add conditional shapes.
        if (length(shapes_present) > 0) {
          p <- p + scale_shape_manual(values = shapes_present)
        }
        if (length(fills_present) > 0) {
          p <- p + scale_fill_manual(values = fills_present)
        }
        if (length(shapes_present) > 0 || length(fills_present) > 0) {
          p <- p + guides(
            shape = guide_legend(title = "Predictions", title.position = "top",
                                 override.aes = list(fill = fills_present, stroke = 1)),
            fill = guide_legend(title = "Predictions", title.position = "top",
                                override.aes = list(shape = shapes_present, size = 2.2))
          )
        }

      } else {
        # Case when not all DBs are in input$db_values_prediction.
        p <- ggplot(prediction_plot, aes(x = C, y = RTsec, color = as.factor(DB),
                                         text = paste0(
                                           "<br><b>", lipidname, "  ---  </b>", adduct_text,

                                           ifelse(IS == "IS","<br><b> - Internal Standard</b>",""),

                                           # Line
                                           "<br>",
                                           paste0(
                                             rep("─", max(1, nchar(
                                               ifelse(IS == "IS",
                                                      " - Internal Standard",
                                                      paste0(ifelse(is.na(lipidname), "", lipidname), "  ", ifelse(is.na(adduct), "", adduct)))
                                             ))),
                                             collapse = ""
                                           ),
                                           "<br><b> - LipidRTid:</b> ", LipidRTid,
                                           "<br><b> - mz:</b> ", round(mz,3),
                                           "<br><b> - RTsec:</b> ", round(RTsec, 3),
                                           "<br><b> - RTsec_rtpred:</b> ", round(RTsec_rtpred, 3),
                                           "<br><b> - Score_rtpred:</b> ", Score_rtpred,
                                           "<br><b> - confidenceLevel:</b> ", confidenceLevel,
                                           "<br><b> - Score_rtpred:</b> ", Score_rtpred,
                                           "<br>   (", checkcandidate, ")",

                                           # other isobars
                                           ifelse(
                                             !is.na(other_isobars_text) & other_isobars_text != "",
                                             paste0(
                                               "<br>",
                                               paste0(
                                                 rep("─", max(1, nchar(
                                                   ifelse(IS == "IS",
                                                          " - Internal Standard",
                                                          paste0(ifelse(is.na(lipidname), "", lipidname),
                                                                 "  ",
                                                                 ifelse(is.na(adduct), "", adduct)))
                                                 ))),
                                                 collapse = ""
                                               ),
                                               other_isobars_text
                                             ),
                                             ""
                                           ),

                                           # highlight
                                           ifelse(
                                             !is.na(comments) & comments != "",
                                             paste0(
                                               "<br>",
                                               paste0(
                                                 rep("─", max(1, nchar(
                                                   ifelse(IS == "IS",
                                                          " - Internal Standard",
                                                          paste0(ifelse(is.na(lipidname), "", lipidname),
                                                                 "  ",
                                                                 ifelse(is.na(adduct), "", adduct)))
                                                 ))),
                                                 collapse = ""
                                               ),
                                               "<br> [ <b>!</b> ] <b>", comments, "</b> [ <b>!</b> ] "
                                             ),
                                             ""
                                           )
                                         ),

                                         alpha = alpha, customdata = index)) +
          geom_point(data = prediction_plot[prediction_plot$shape == 16,], size = 2.2, stroke = 0.6) +
          geom_point(data = prediction_plot[prediction_plot$shape == 22,],
                     aes(shape = legend_label, fill = legend_label), size = 2.2, stroke = 0.6) +

          labs(x = "Carbon number", y = "Retention time (s)") +
          scale_color_manual(name = "Double \nbonds",
                             values = c("#212121", "#EE2C2C", "#FF9500", "#FFE100",
                                        "#90EE90", "#228B22", "#B0E2FF", "#00B2EE",
                                        "#3A5FCD", "#BF3EFF", "#68228B", "#795548",
                                        "#9E9E9E", "#FFB6C1", "#FFD39B", "#FFF68F",
                                        "#B8860B", "#8B864E", "#B87F9E", "#7FFF00",
                                        "#FF69B4", "#DEBF4B", "#51FFFC", "#EAC4FF",
                                        "#D9FF68",
                                        "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
                                        "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
                                        "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173",
                                        "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
                                        "#969696", "#DD1C77", "#6BAED6", "#74C476", "#9E9AC8",
                                        "#E7969C", "#CEDB9C", "#9C9EDE", "#DBDB8D", "#C49C94",
                                        "#F7B6D2", "#C7C7C7", "#AEC7E8", "#FFBB78", "#98DF8A",
                                        "#FF9896", "#C5B0D5", "#C49C94", "#FDD0A2", "#BDBDBD",
                                        "#5254A3", "#6B6ECF", "#9C9EDE", "#637939", "#B5CF6B",
                                        "#CEDB9C", "#8CA252", "#BD9E39", "#E7BA52", "#E7CB94",
                                        "#843C39", "#AD494A", "#D6616B", "#E7969C", "#7B4173",
                                        "#A55194", "#CE6DBD", "#DE9ED6", "#3182BD", "#6BAED6",
                                        "#9ECAE1", "#C6DBEF", "#E6550D", "#FD8D3C", "#FDAE6B",
                                        "#FDD0A2", "#31A354", "#74C476", "#A1D99B", "#C7E9C0",
                                        "#756BB1", "#9E9AC8", "#BCBDDC", "#DADAEB", "#636363",
                                        "#969696", "#BDBDBD", "#D9D9D9", "#393B79", "#5254A3",
                                        "#6B6ECF", "#9C9EDE", "#637939", "#8CA252", "#B5CF6B",
                                        "#CEDB9C", "#8C6D31", "#BD9E39", "#E7BA52", "#E7CB94"),
                             breaks = seq(0:24)-1) +
          scale_shape_manual(values = shapes_present) +
          scale_fill_manual(values = fills_present) +
          guides(shape = guide_legend(title = "Predictions", title.position = "top",
                                      override.aes = list(fill = fills_present, stroke = 1)),
                 fill = guide_legend(title = "Predictions", title.position = "top",
                                     override.aes = list(shape = shapes_present, size = 2.2))) +
          scale_x_continuous(limits = c(minc_plot, maxc_plot), breaks = seq(minc_plot, maxc_plot, by = 2)) +
          scale_y_continuous(limits = c(minrt_plot, maxrt_plot)) +
          theme_minimal() +
          theme(legend.position = "none") +

          geom_point(data = subset(prediction_plot, !is.na(comments) & comments != ""),
                     aes(x = C, y = RTsec), shape = 22, size = 5,
                     fill = "black", alpha = 0.25, color = NA, inherit.aes = FALSE)
      }

      # Add subtitle with lipid classes
      # Include which classes were used to build the models ....................
      models <- rtdata$prediction$models
      predict_info <- rtdata$prediction$information
      lipid_classes <- NULL
      if (!is.null(models[[lipidclass]]$model$data) && !is.null(models[[lipidclass]]$model$data$LipidSubclass)) {
        lipid_classes <- unique(models[[lipidclass]]$model$data$LipidSubclass)
      } else if (!is.null(models[[lipidclass]]$model$model) && !is.null(models[[lipidclass]]$model$model$LipidSubclass)) {
        lipid_classes <- unique(models[[lipidclass]]$model$model$LipidSubclass)
      } else {
      }

      subtitle_text <- character(0)
      if (!is.null(lipid_classes)) {
        if (length(lipid_classes) == 1) {
          subtitle_text <- paste0("This model was built using the lipid subclass: ", lipid_classes, ".")
        } else {
          subtitle_text <- paste0("This model was built using the following lipid subclasses: ",
                                 paste0(lipid_classes, collapse = ", "), ".")
        }

        if (lipidclass %in% unique(rtdata$prediction$candidates$LipidSubclass)) {
          subtitle_text <- paste(
            subtitle_text, "\n ·",
            ifelse(is.na(predict_info$new_rtannotations[predict_info$LipidSubclass == lipidclass]),
                   0,
                   predict_info$new_rtannotations[predict_info$LipidSubclass == lipidclass]),
            "accepted candidates,  ",
            ifelse(is.na(predict_info$discarded_rtannotations[predict_info$LipidSubclass == lipidclass]),
                   0,
                   predict_info$discarded_rtannotations[predict_info$LipidSubclass == lipidclass]),
            "discarded candidates.")
        } else {
          subtitle_text <- paste(
            subtitle_text, "\n ·","No proposed candidates for this lipid class."
          )
        }


      } else {
        subtitle_text <- "No model for this lipid subclass."
        subtitle_text <- paste(
          subtitle_text, "\n ·","0 accepted candidates,  ", "0 discarded candidates."
        )
      }



      p <- p +
        ggtitle(subtitle_text) +
        theme(plot.title = element_text(size = 9, color = "grey41"))

      # If there is a prediction_layer, plot it ................................
      if (exists("prediction_layer")) {
        if (length(prediction_layer) != 0) {
          p <- p +
            prediction_layer +
            theme(legend.position = "none")
        }
      }

      # Render with Plotly .....................................................
      plotly::layout(
        plotly::ggplotly(
          p,
          tooltip = c("text"),
          source = "manualPredictionPlot_source"
        ),
        dragmode = "lasso"
      )
    }
  }


  #============================================================================#
  # Render the initial plot
  #============================================================================#
  output$manualPredictionPlot <- plotly::renderPlotly({
    rtdata <- rtdataReactive()
    if (!is.null(rtdata$prediction$candidates)) {
      rtdata$prediction$candidates <-
        create_other_isobars(rtdata$prediction$candidates,
                             typeanalysis = "prediction")
    }
    rtdataReactive(rtdata)

    renderPredictionGraph()
  })


  #============================================================================#
  # Downloads
  #============================================================================#
  # Download BATCH RESULTS .....................................................
  output$downloadPeaklistResults <- downloadHandler(
    filename = function(){paste(input$jobname, "BatchResults.zip", sep="_")},
    content = function(file){

        notif_id <- showNotification(HTML("Downloading 'Peaklist Results'... <span class='loader'></span>"),
                                     type = "default", duration = NULL)

        rtdata <- rtdataReactive()

        # go to a temp dir to avoid permission issues
        owd <- setwd(tempdir())
        on.exit(setwd(owd));
        files <- NULL

        # extract data
        peaklist <- rtdata$features
        summarize <- rtdata$summarizedFeatures

        # write peaklist files
        write.csv(peaklist, "Peaklist.csv", row.names = FALSE)
        files <- c(files, "Peaklist.csv")

        # write summarize files
        write.csv(summarize, "summarizePeaklist.csv", row.names = FALSE)
        files <- c(files, "summarizePeaklist.csv")

        # create the zip file
        zip::zip(file, files)

        removeNotification(notif_id)
        showNotification("'Peaklist Results' has been successfully downloaded.", type = "message", duration = NULL)
    }
  )

  # DOWNLOAD FILTER PLOTS .....................................................
  output$downloadFilterResults <- downloadHandler(
    filename = function() {
      paste(input$jobname, "FilterResults.zip", sep = "_")
    },

    content = function(file) {

      notif_id <- showNotification(
        HTML("Downloading 'Filter Results'... <span class='loader'></span>"),
        type = "default", duration = NULL
      )

      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      files <- NULL

      rtdata <- rtdataReactive()

      ## -------- PDF with filtering plots --------
      fileName_pdf <- paste(input$jobname, "FilterPlots.pdf", sep = "_")
      grDevices::pdf(file = fileName_pdf, width = 20, height = 12)

      for (p in seq_along(rtdata$filtering$plots)) {
        print(rtdata$filtering$plots[[p]])
      }

      grDevices::dev.off()
      files <- c(files, fileName_pdf)

      ## -------- TXT with model report --------
      model_report <- LipidRT::modelreport(rtdata, "filtering")
      fileName_report <- paste(input$jobname, "FilterReportModel.txt", sep = "_")
      writeLines(model_report, con = fileName_report)
      files <- c(files, fileName_report)

      ## -------- TXT with filtering summary --------
      fileName_filteringInfo <- paste(
        input$jobname, "FilteringSummary.txt", sep = "_"
      )

      writeLines(LipidRT::get_filtering_summary(rtdata), con = fileName_filteringInfo)

      files <- c(files, fileName_filteringInfo)

      ## -------- ZIP all files --------
      if (length(files) > 0) {
        zip::zip(zipfile = file, files = files)

        removeNotification(notif_id)
        showNotification(
          "'Filter Results' has been successfully downloaded.",
          type = "message", duration = NULL
        )
      } else {
        removeNotification(notif_id)
        showNotification(
          "No files were generated for download.",
          type = "error", duration = NULL
        )
      }
    }
  )


  # Download PREDICTION PLOTS .................................................
  output$downloadPredictionResults <- downloadHandler(

    filename = function() {
      paste(input$jobname, "PredictionResults.zip", sep = "_")
    },

    content = function(file) {

      notif_id <- showNotification(
        HTML("Downloading 'Prediction Results'... <span class='loader'></span>"),
        type = "default", duration = NULL
      )

      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      files <- NULL

      rtdata <- rtdataReactive()

      ## -------- PDF with prediction plots --------
      fileName_pdf <- paste(input$jobname, "PredictionPlots.pdf", sep = "_")
      grDevices::pdf(file = fileName_pdf, width = 20, height = 12)

      for (p in seq_along(rtdata$prediction$plots)) {
        print(rtdata$prediction$plots[[p]])
      }

      grDevices::dev.off()
      files <- c(files, fileName_pdf)

      ## -------- TXT model report --------
      model_report <- LipidRT::modelreport(rtdata, "prediction")
      fileName_report <- paste(input$jobname, "PredictionReportModel.txt", sep = "_")
      writeLines(model_report, con = fileName_report)
      files <- c(files, fileName_report)

      ## -------- TXT prediction summary --------
      fileName_predictionInfo <- paste(
        input$jobname, "PredictionSummary.txt", sep = "_"
      )

      writeLines(LipidRT::get_prediction_summary(rtdata), con = fileName_predictionInfo)
      files <- c(files, fileName_predictionInfo)

      ## -------- ZIP all files --------
      if (length(files) > 0) {
        zip::zip(zipfile = file, files = files)

        removeNotification(notif_id)
        showNotification(
          "'Prediction Results' has been successfully downloaded.",
          type = "message", duration = NULL
        )
      } else {
        removeNotification(notif_id)
        showNotification(
          "No files were generated for download.",
          type = "error", duration = NULL
        )
      }
    }
  )


  # Download TARGET RESULTS ....................................................
  output$downloadTargetResults <- downloadHandler(
    filename = function(){paste(input$jobname, "TargetPrediction.zip", sep="_")},
    content = function(file){

      notif_id <- showNotification(HTML("Downloading 'Target Predictions'... <span class='loader'></span>"),
                                   type = "default", duration = NULL)

      # go to a temp dir to avoid permission issues
      owd <- setwd(tempdir())
      on.exit(setwd(owd));

      # extract data
      target <- targetpredictionsReactive()

      # write files
      write.csv(target, "TargetPrediction.csv", row.names = FALSE)

      files <- c("TargetPrediction.csv")

      # create the zip file
      zip::zip(file, files)

      removeNotification(notif_id)
      showNotification("'Target Results' has been successfully downloaded.", type = "message", duration = NULL)
    }
  )



  # Download .RDATA ...........................................................
  output$download_rtdata1 <- downloadHandler(
    filename = function() {
      rtdata <- rtdataReactive()

      if (is.null(rtdata)) {
        return(NULL)  # If there is no data, skip the download.
      } else {
        paste(input$jobname, "rtdata.RData", sep = "_")
      }
    },
    content = function(file) {

      notif_id <- showNotification(HTML("Downloading 'rtdata' object... <span class='loader'></span>"),
                                   type = "default", duration = NULL)

      rtdata <- rtdataReactive()

      if (is.null(rtdata)) {
        removeNotification(notif_id)
        howNotification("No rtdata object were generated for download.", type = "error", duration = NULL)
        return()  # If there is no data, do nothing.
      } else {
          rtdata <- rtdata

          if (!is.null(rtdata$prediction$candidates$other_isobars)) {
            rtdata$prediction$candidates$other_isobars <- NULL
          }

          save(rtdata, file = file)
          removeNotification(notif_id)
          showNotification("'rtdata' object has been successfully downloaded.", type = "message", duration = NULL)
      }
    }
  )

  output$download_rtdata2 <- downloadHandler(
    filename = function() {
      rtdata <- rtdataReactive()

      if (is.null(rtdata)) {
        return(NULL)  # If there is no data, skip the download.
      } else {
        paste(input$jobname, "rtdata.RData", sep = "_")
      }
    },
    content = function(file) {

      notif_id <- showNotification(HTML("Downloading 'rtdata' object... <span class='loader'></span>"),
                                   type = "default", duration = NULL)

      rtdata <- rtdataReactive()

      if (is.null(rtdata)) {
        removeNotification(notif_id)
        howNotification("No rtdata object were generated for download.", type = "error", duration = NULL)
        return()  # If there is no data, do nothing.
      } else {
        rtdata <- rtdata

        if (!is.null(rtdata$prediction$candidates$other_isobars)) {
          rtdata$prediction$candidates$other_isobars <- NULL
        }

        save(rtdata, file = file)
        removeNotification(notif_id)
        showNotification("'rtdata' object has been successfully downloaded.", type = "message", duration = NULL)
      }
    }
  )


  ##############################################################################
  # Extra Functions
  ##############################################################################
  # select_lipidclasses <- function(analysis_step, type_lipidclasses) {
  #   polarity <- rtdataReactive()$metadata$polarity
  #
  #   adductsTable <- adductsTableReactive()
  #   if (!is.null(polarity)) {
  #     if (polarity == "positive") {
  #       valid_adducts <- adductsTable$adduct[adductsTable$charge > 0]
  #     } else if (polarity == "negative") {
  #       valid_adducts <- adductsTable$adduct[adductsTable$charge < 0]
  #     } else if (polarity == "both") {
  #       valid_adducts <- c()
  #     }
  #   } else {
  #     valid_adducts <- adductsTable$adduct
  #   }
  #
  #   dbs <- dbsReactive()
  #
  #   # Base labels according to the analysis step
  #   label <- switch(
  #     analysis_step,
  #     "lipidclasses_predict" = "Lipid classes to predict",
  #     "lipidclasses_annotate" = "Lipid classes to annotate"
  #   )
  #
  #   label_custom <- paste(label, "(empty list)")
  #   label_nopolarity <- "No data imported"
  #   label_nofile <- "No file selected"
  #
  #   # Determine if there is data or a file
  #   no_polarity <- is.null(polarity)
  #   both_polarity <- polarity == "both"
  #   no_file <- is.null(dbs)
  #
  #   # If there is no polarity and it is the default list, show a "no data" message
  #   if (no_polarity && type_lipidclasses == "lipidclasses_default") {
  #     label_html <- HTML(paste0(
  #       "<b>", label_custom, "</b><br>",
  #       "<span style='color:gray;'>&nbsp;&nbsp;&nbsp;&nbsp;", label_nopolarity, "</span>"
  #     ))
  #
  #     updateCheckboxGroupInput(session, analysis_step,
  #                              label = label_html,
  #                              choices = list(),
  #                              selected = NULL)
  #
  #     return()
  #   }
  #
  #   if (both_polarity) {
  #     label_html <- HTML(paste0(
  #       "<b>", label_custom, "</b><br>",
  #       "<span style='color:gray;'>&nbsp;&nbsp;&nbsp;&nbsp;
  #       Adducts must correspond to either positive or negative polarity only.\n
  #       Please ensure all adducts belong to a single polarity group.</span>"
  #     ))
  #
  #     updateCheckboxGroupInput(session, analysis_step,
  #                              label = label_html,
  #                              choices = list(),
  #                              selected = NULL)
  #     return()
  #   }
  #
  #   # Build options if there is a database
  #   if (!no_file) {
  #
  #     choicesNames <- list()
  #     choicesValues <- list()
  #
  #     idx <- 1
  #
  #     for (i in seq_along(dbs)) {
  #
  #       lipid <- names(dbs)[i]
  #
  #       if (lipid != "adductsTable") {
  #
  #         adducts_valid <- intersect(
  #           colnames(dbs[[lipid]])[-c(1:3)],
  #           valid_adducts
  #         )
  #
  #         # Only if there are valid adducts
  #         if (length(adducts_valid) > 0) {
  #
  #           choicesNames[[idx]] <- HTML(paste0(
  #             "<b style='color:black;'>", lipid, "&nbsp;&nbsp;&nbsp;&nbsp;</b>",
  #             "<small>",
  #             "<span style='color:gray;'>( </span>",
  #             paste(
  #               "<span style='color:#669999;'>",
  #               adducts_valid,
  #               "</span>",
  #               collapse = paste0(
  #                 "<span style='color:gray;'>&nbsp;",
  #                 gsub("\\\\", "", "\\|"),
  #                 "&nbsp;</span>"
  #               )
  #             ),
  #             "<span style='color:gray;'> )</span>",
  #             "</small>"
  #           ))
  #
  #           choicesValues[[idx]] <- lipid
  #           idx <- idx + 1
  #         }
  #       }
  #     }
  #
  #     lipid_adducts <- choicesNames
  #     names(lipid_adducts) <- unlist(choicesValues)
  #
  #     choicesNamesReactive(lipid_adducts)
  #   }
  #
  #   # Add polarity information if available
  #   if (!no_polarity && type_lipidclasses == "lipidclasses_default") {
  #     if (polarity == "positive") {
  #       label <- paste(label, "for ESI + :")
  #     } else if (polarity == "negative") {
  #       label <- paste(label, "for ESI - :")
  #     }
  #   }
  #
  #   # Show options according to the list type
  #   if (type_lipidclasses == "lipidclasses_default") {
  #     label <- gsub("\\(empty list\\)", "", label)
  #
  #     if (!no_polarity && !no_file) {
  #       label <- paste0(gsub(":", "", label), " (Default list) :")
  #       updateCheckboxGroupInput(session, analysis_step,
  #                                label = label,
  #                                choiceNames = choicesNames,
  #                                choiceValues = choicesValues,
  #                                selected = choicesValues)
  #     } else {
  #       label_html <- HTML(paste0(
  #         "<b>", label_custom, "</b><br>",
  #         "<span style='color:gray;'>&nbsp;&nbsp;&nbsp;&nbsp;",
  #         if (no_polarity) label_nopolarity else label_nofile,
  #         "</span>"
  #       ))
  #
  #       updateCheckboxGroupInput(session, analysis_step,
  #                                label = label_html,
  #                                choices = list(),
  #                                selected = NULL)
  #     }
  #
  #   } else if (type_lipidclasses == "lipidclasses_custom") {
  #     label <- gsub("\\(empty list\\)", "", label)
  #
  #     if (!no_file) {
  #       label <- paste0(gsub(":", "", label), " (Custom list) :")
  #       updateCheckboxGroupInput(session, analysis_step,
  #                                label = label,
  #                                choiceNames = choicesNames,
  #                                choiceValues = choicesValues,
  #                                selected = choicesValues)
  #     } else {
  #       label_html <- HTML(paste0(
  #         "<b>", label_custom, "</b><br>",
  #         "<span style='color:gray;'>&nbsp;&nbsp;&nbsp;&nbsp;", label_nofile, "</span>"
  #       ))
  #
  #       updateCheckboxGroupInput(session, analysis_step,
  #                                label = label_html,
  #                                choices = list(),
  #                                selected = NULL)
  #     }
  #   }
  # }
  select_lipidclasses <- function(analysis_step, type_lipidclasses) {

    rtdata <- rtdataReactive()
    dbs <- dbsReactive()

    # =========================
    # EARLY EXIT
    # =========================
    if (is.null(rtdata) || is.null(dbs) || is.null(rtdata$metadata$polarity)) {

      shinyWidgets::updatePickerInput(
        session,
        analysis_step,
        choices = list(),
        selected = NULL
      )

      return()
    }

    polarity <- rtdata$metadata$polarity
    adductsTable <- adductsTableReactive()

    # =========================
    # VALID ADDUCT FILTER
    # =========================
    if (!is.null(polarity)) {
      if (polarity == "positive") {
        valid_adducts <- adductsTable$adduct[adductsTable$charge > 0]
      } else if (polarity == "negative") {
        valid_adducts <- adductsTable$adduct[adductsTable$charge < 0]
      } else {
        valid_adducts <- c()
      }
    } else {
      valid_adducts <- adductsTable$adduct
    }

    no_file <- is.null(dbs)

    # =========================
    # BUILD SIMPLE LIST (ONLY LIPIDS)
    # =========================
    lipid_choices <- c()

    if (!no_file) {

      for (i in seq_along(dbs)) {

        lipid <- names(dbs)[i]

        if (lipid != "adductsTable") {

          adducts_valid <- intersect(
            colnames(dbs[[lipid]])[-c(1:3)],
            valid_adducts
          )

          # only include lipid if it has valid adducts
          if (length(adducts_valid) > 0) {
            lipid_choices <- c(lipid_choices, lipid)
          }
        }
      }
    }

    # named vector (pickerInput format)
    lipid_choices <- setNames(lipid_choices, lipid_choices)

    # =========================
    # UPDATE PICKER
    # =========================
    if (!no_file) {

      shinyWidgets::updatePickerInput(
        session,
        analysis_step,
        choices = lipid_choices,
        selected = lipid_choices
      )

    } else {

      shinyWidgets::updatePickerInput(
        session,
        analysis_step,
        choices = list(),
        selected = NULL
      )
    }
  }


  # Function to extract the databases (DB) associated with each lipid class.
  # This is used for plotting and also for generating the LipidClasses list in
  # annotation and prediction.
  DBextraction <- function(rtdata) {

    # Helper function to get sorted unique values grouped by LipidClass.
    get_unique_sorted <- function(data, key) {
      split_values <- split(data[[key]], data$LipidSubclass)
      lapply(split_values, function(x) sort(unique(x)))
    }

    dbChoices <- list()  # Initialize the list to store DBs by context (filtering, prediction, none)

    filtering <- FALSE
    prediction <- FALSE

    # Check for filtering and prediction availability
    if (!is.null(rtdata)) {
      if (!is.null(rtdata$filtering$data)) filtering <- TRUE
      if (!is.null(rtdata$prediction$candidates)) prediction <- TRUE
    }

    if (filtering | prediction) {

      # Handle filtering data
      if (filtering) {
        valid_filter_data <- subset(rtdata$filtering$data,
                                    !is.na(LipidSubclass) &
                                      LipidSubclass != "" &
                                      !(checkfilter %in% c("Ignore duplicated", "Not pass QC",
                                                      "Unsupported class"))&
                                      !is.na(checkfilter))

        db_initial_filter <- lapply(
          split(valid_filter_data$DB, valid_filter_data$LipidSubclass),
          function(x) sort(unique(x))
        )

        if (!is.list(dbChoices[["filtering"]])) dbChoices[["filtering"]] <- list()
        dbChoices[["filtering"]] <- db_initial_filter
      }

      # Handle prediction data
      if (prediction) {
        if (nrow(rtdata$prediction$candidates) != 0) {

        candidates <- subset(rtdata$prediction$candidates,
                             !is.na(LipidSubclass) &
                               LipidSubclass != "" &
                               !(checkcandidate %in% c("Removed")) &
                               !is.na(checkcandidate))
        filtered <- subset(rtdata$filtering$data,
                           !is.na(LipidSubclass) &
                             LipidSubclass != "" &
                             !(checkfilter %in% c("Ignore duplicated", "Not pass QC",
                                             "Unsupported class")) &
                             !is.na(checkfilter))

        db_candidates <- get_unique_sorted(candidates, "DB")
        db_filtered <- get_unique_sorted(filtered, "DB")

        lipid_classes_prediction <- names(rtdata$prediction$plots)

        db_initial_prediction <- Map(function(candidates, filtered) {
          sort(unique(c(candidates, filtered)))
        }, db_candidates[lipid_classes_prediction], db_filtered[lipid_classes_prediction])

        names(db_initial_prediction) <- lipid_classes_prediction

        if (!is.list(dbChoices[["prediction"]])) dbChoices[["prediction"]] <- list()
        dbChoices[["prediction"]] <- db_initial_prediction
        }
      } else {
        if (filtering) {
          valid_filter_data <- subset(rtdata$filtering$data,
                                      !is.na(LipidSubclass) &
                                        LipidSubclass != "" &
                                        !(checkfilter %in% c("Ignore duplicated", "Not pass QC",
                                                             "Unsupported class")))

          db_initial_filter <- lapply(
            split(valid_filter_data$DB, valid_filter_data$LipidSubclass),
            function(x) sort(unique(x))
          )

          if (!is.list(dbChoices[["filtering"]])) dbChoices[["filtering"]] <- list()
          dbChoices[["filtering"]] <- db_initial_filter
        }
      }

    } else {
      separate_features <- subset(rtdata$data, !is.na(LipidSubclass) & LipidSubclass != "")

      db_values <- lapply(
        split(separate_features$DB, separate_features$LipidSubclass),
        function(x) sort(unique(x))
      )

      if (!is.list(dbChoices[["none"]])) dbChoices[["none"]] <- list()
      dbChoices[["none"]] <- db_values
    }

    return(dbChoices)  # Return the list of DBs grouped by context and lipid class
  }



  create_other_isobars <- function(df, typeanalysis) {

    Filter_final <- FALSE
    if ("Filter_final" %in% colnames(df)) {
      Filter_final <- TRUE
      required_cols <- c("LipidRTid", "lipidname", "adduct", "RTsec", "RTsec_rtpred",
                         "Score_rtpred", "Filter_final")
    } else {
      if (typeanalysis == "prediction") {
        required_cols <- c("LipidRTid", "lipidname", "adduct", "RTsec", "RTsec_rtpred",
                           "Score_rtpred", "checkcandidate")
      } else if (typeanalysis == "filtering") {
        required_cols <- c("LipidRTid", "lipidname", "adduct", "RTsec", "RTsec_rtpred",
                           "Score_rtpred", "checkfilter")
      }
    }


      if (!all(required_cols %in% names(df))) {
        stop("Missing required columns to create other_isobars")
      }

      if ("other_isobars" %in% names(df)) {
        df$other_isobars <- NULL
      }

      # Ensure character type (only for processing, does not change logic)
      df$lipidname      <- as.character(df$lipidname)
      df$adduct         <- as.character(df$adduct)
      df$Score_rtpred   <- as.character(df$Score_rtpred)

      # --------------------------------------------------
      # 1) Use ONLY Kept / Discarded as source
      # --------------------------------------------------
      # if (typeanalysis == "prediction") {
      #   df$checkcandidate <- as.character(df$checkcandidate)
      #   source_df <- df[df$checkcandidate %in% c("Kept", "Discarded"), , drop = FALSE]
      # } else if (typeanalysis == "filtering") {
      #   df$checkfilter <- as.character(df$checkfilter)
      #   source_df <- df[!df$checkfilter %in% c("Ignore duplicated"), , drop = FALSE]
      # }
      if (typeanalysis == "prediction") {
        if (Filter_final) {
          df$Filter_final <- as.character(df$Filter_final)
          source_df <- df[df$Filter_final %in% c("Kept", "Kept manually by the user", "Recovered with RT modelling"), , drop = FALSE]
        } else {
          df$checkcandidate <- as.character(df$checkcandidate)
          source_df <- df[df$checkcandidate %in% c("Kept", "Kept manually by the user", "Recovered with RT modelling"), , drop = FALSE]
        }
      } else if (typeanalysis == "filtering") {
        if (Filter_final) {
          df$Filter_final <- as.character(df$Filter_final)
          source_df <- df[df$Filter_final %in% c("Kept", "Kept manually by the user", "Recovered with RT modelling"), , drop = FALSE]
        } else {
          df$checkfilter <- as.character(df$checkfilter)
          source_df <- df[df$checkfilter %in% c("Kept", "Kept manually by the user", "Recovered with RT modelling"), , drop = FALSE]
        }
      }



      # If no valid rows exist, create empty column and return
      if (nrow(source_df) == 0) {
        df$other_isobars <- rep("", nrow(df))
        return(df)

      }

      # Isobar identifier used only in the source data
      source_df$.isobar_key <- paste(
        source_df$lipidname,
        source_df$adduct,
        source_df$RTsec,
        source_df$RTsec_rtpred,
        source_df$Score_rtpred,
        sep = ","
      )

      # Group isobars by LipidRTid
      iso_by_rt <- split(source_df$.isobar_key, source_df$LipidRTid)

      # --------------------------------------------------
      # 2) Assign other_isobars to ALL rows
      # --------------------------------------------------
      df$other_isobars <- mapply(
        function(rt, lipid, adduct, rtsec, rtpred, score) {

          key <- paste(
            lipid,
            adduct,
            rtsec,
            rtpred,
            score,
            sep = ","
          )

          if (!as.character(rt) %in% names(iso_by_rt)) {
            return("")
          }

          others <- setdiff(unique(iso_by_rt[[as.character(rt)]]), key)

          if (length(others) == 0) {
            ""
          } else {
            paste(others, collapse = "|")
          }
        },
        df$LipidRTid,
        df$lipidname,
        df$adduct,
        df$RTsec,
        df$RTsec_rtpred,
        df$Score_rtpred,
        USE.NAMES = FALSE
      )

    return(df)
  }



  session$onSessionEnded(function() {
    stopApp()
  })


}

# Run the application
shinyApp(ui = ui, server = server)

