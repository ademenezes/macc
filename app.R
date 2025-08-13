library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(janitor)
library(tidyr)
library(mc2d)
library(scales)
library(gridExtra)
library(plotly)


   safe_q <- function(x, p) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    as.numeric(stats::quantile(x, p, type = 7))
  }

# Extended Monte Carlo-enhanced MACC computation with adoption and energy uncertainty
  compute_costs_macc_uncertain_extended <- function(
    sel, r, subsidy_pct = 0, tariff_price = 0,
    adoption_rate_mean = 0.3, adoption_rate_sd = 0.05,
    energy_saving_sd_frac = 0.2, n_sim = 1000,
    benefits = list(ws=TRUE, es=TRUE, yi=TRUE)
  ) {
    s      <- subsidy_pct / 100
    r_eff  <- if (isTRUE(is.na(r)) || r <= 0) 1e-9 else r  # avoid div-by-zero
    # ensure timeframe positive in formulas
    sel$timeframe <- ifelse(is.na(sel$timeframe) | sel$timeframe <= 0, 1, sel$timeframe)

    sims <- lapply(seq_len(n_sim), function(i) {
      adoption_rate_sim <- rnorm(1, mean = adoption_rate_mean, sd = adoption_rate_sd)

      sel %>%
        mutate(
          alpha       = alpha_pct / 100,
          WS          = ws_base_km3 * 1000,

          c_ha_sim    = rnorm(n(), mean = coalesce(c_ha, 0), sd = 0.1 * pmax(coalesce(c_ha, 0), 1e-9)),
          b_ws_sim    = rnorm(n(), mean = coalesce(b_ws, 0), sd = 0.2 * pmax(coalesce(b_ws, 0), 1e-9)),
          b_es_sim    = rnorm(n(), mean = coalesce(b_es, 0), sd = 0.2 * pmax(coalesce(b_es, 0), 1e-9)),
          b_yi_sim    = rnorm(n(), mean = coalesce(b_yi, 0), sd = 0.2 * pmax(coalesce(b_yi, 0), 1e-9)),
          es_mwh_sim  = rnorm(n(), mean = coalesce(es_mwh, 0), sd = energy_saving_sd_frac * pmax(coalesce(es_mwh, 0), 1e-9)),

          alpha_sim   = pmin(pmax(alpha * (adoption_rate_sim / adoption_rate_mean), 0), 1),

          C_tot       = coalesce(c_tot_in, a_pot * c_ha_sim * alpha_sim),
          C_sub       = (1 - s) * C_tot,

          b_ws_eff    = if (isTRUE(benefits$ws)) b_ws_sim else 0,
          b_es_eff    = if (isTRUE(benefits$es)) b_es_sim else 0,
          b_yi_eff    = if (isTRUE(benefits$yi)) b_yi_sim else 0,
          b_total_eff = b_ws_eff + b_es_eff + b_yi_eff,

          PV_benefit  = (b_total_eff / r_eff) * (1 - (1 + r_eff)^(-timeframe)),
          Net_Cost    = C_sub - PV_benefit,
          EAC         = (r_eff * Net_Cost) / (1 - (1 + r_eff)^(-timeframe)),

          # Avoid division by zero → set MC to NA when WS <= 0
          MC          = ifelse(WS > 0, EAC / WS, NA_real_),

          scenario    = i
        ) %>%
        select(measure, WS, Net_Cost, MC, scenario)
    })

    bind_rows(sims) %>%
      group_by(measure) %>%
      summarise(
        WS        = safe_q(WS, 0.5),        # median WS (or NA if none)
        Net_Cost  = safe_q(Net_Cost, 0.5),  # median Net Cost
        MC_median = safe_q(MC, 0.5),
        MC_p05    = safe_q(MC, 0.05),
        MC_p95    = safe_q(MC, 0.95),
        .groups = "drop"
      ) %>%
      arrange(MC_median)
  }



# Settings
PARAM_FILE <- "params_3.csv"
dt         <- 1

# Read and clean
raw <- read.csv(PARAM_FILE, stringsAsFactors = FALSE)
df  <- raw %>% clean_names()

# Mapping
mapping_clean <- list(
  country     = "^region$",
  level       = "level_of_applicability|level",
  measure     = "^measure$",
  a_pot       = "potential_area",
  ews_pct     = "expected_water_saving",
  wr          = "water_requirement",
  c_ha        = "cost_per_hectare",
  alpha_pct   = "adaptation_coefficient",
  ws_base_km3 = "water_saving.*1000",
  es_mwh      = "energy_saving.*m.?wh",
  yi_pct      = "^yield_increase$",
  b_ws        = "water_saving_benefit",
  b_es        = "energy_saving_benefit",
  b_yi        = "yield_increase_benefit",
  b_total     = "total_benefit",
  t_ws        = "government_target",
  c_tot_in    = "tentative_total_cost",
  timeframe   = "^timeframe$",
  bank        = "bank_involvement",
  water_price = "water_price"
)
for(short in names(mapping_clean)) {
  pat   <- mapping_clean[[short]]
  found <- grep(pat, names(df), ignore.case = TRUE, value = TRUE)
  if(length(found) == 1) {
    df <- df %>% rename(!!short := !!sym(found))
  } else {
    warning(sprintf("Mapping for '%s' found %d: %s",
                    short, length(found), paste(found, collapse = ", ")))
  }
}

# Convert numeric
num_vars <- c("a_pot","ews_pct","wr","c_ha","alpha_pct",
              "ws_base_km3","es_mwh","yi_pct",
              "b_ws","b_es","b_yi","b_total","t_ws",
              "c_tot_in","timeframe","water_price")
present_nums <- intersect(num_vars, names(df))
df <- df %>% mutate(across(all_of(present_nums), ~ as.numeric(gsub("[^0-9\\.]", "", .))))

# MACC function
compute_costs_macc <- function(sel, r, subsidy_pct = 0, tariff_price = 0, benefits = list(ws=TRUE, es=TRUE, yi=TRUE)) {
  horizon <- sel$timeframe[1]
  s       <- subsidy_pct / 100

  sel %>%
    mutate(
      alpha = alpha_pct / 100,
      WS    = ws_base_km3 * 1000,

      # Total (upfront) cost
      C_tot = case_when(
        !is.na(c_tot_in) & c_tot_in > 0 ~ c_tot_in,
        TRUE                            ~ a_pot * c_ha * alpha
      ),
      C_sub = C_tot * (1 - s),

      # Build total benefit from selected components only
      b_ws_eff = if (isTRUE(benefits$ws)) coalesce(b_ws, 0) else 0,
      b_es_eff = if (isTRUE(benefits$es)) coalesce(b_es, 0) else 0,
      b_yi_eff = if (isTRUE(benefits$yi)) coalesce(b_yi, 0) else 0,
      b_total_eff = b_ws_eff + b_es_eff + b_yi_eff,

      PV_benefit = (b_total_eff / r) * (1 - (1 + r)^(-horizon)),
      Net_Cost   = C_sub - PV_benefit,
      EAC        = (r * Net_Cost) / (1 - (1 + r)^(-horizon)),
      MC         = EAC / WS,

      P_eff     = pmax(water_price, tariff_price),
      price_gap = MC - P_eff,
      is_cost_effective = MC <= P_eff
    ) %>%
    arrange(MC) %>%
    mutate(cumWS = cumsum(WS))
}


simulate_area <- function(A_pot, adoption_rate, midpoint_year, horizon, dt = 1) {
  years <- seq(0, horizon, by = dt)
  A <- A_pot / (1 + exp(-adoption_rate * (years - midpoint_year)))
  data.frame(year = years, area = A)
}

simulate_energy_saving <- function(area_df, es_mwh_total, A_pot) {
  area_df %>% mutate(energy_saved = (area / A_pot) * es_mwh_total)
}

simulate_yield_benefit <- function(area_df, B_yi_total, A_pot) {
  area_df %>% mutate(yield_benefit = (area / A_pot) * B_yi_total)
}

simulate_total_benefit <- function(area_df, B_ws, B_es, B_yi, A_pot, benefits = list(ws=TRUE, es=TRUE, yi=TRUE)) {
  B_ws_eff <- if (isTRUE(benefits$ws)) B_ws else 0
  B_es_eff <- if (isTRUE(benefits$es)) B_es else 0
  B_yi_eff <- if (isTRUE(benefits$yi)) B_yi else 0

  area_df %>%
    mutate(
      b_ws_t = (area / A_pot) * B_ws_eff,
      b_es_t = (area / A_pot) * B_es_eff,
      b_yi_t = (area / A_pot) * B_yi_eff,
      b_total = b_ws_t + b_es_t + b_yi_t
    )
}


# UI
ui <- fluidPage(
  
  # --- Top-right logo ---
  tags$head(
    tags$style(HTML("
      #logo-container {
        position: fixed;        /* stays in the top-right even when scrolling */
        top: 10px;
        right: 10px;
        z-index: 1000;          /* above other UI */
      }
      #logo-container img {
        height: 50px;           /* tweak as you like */
      }
      @media (max-width: 768px) {
        #logo-container img { height: 40px; } /* smaller on mobile */
      }
    "))
  ),
  tags$div(
    id = "logo-container",
    # Make it clickable? wrap in tags$a(href="https://your.link", target="_blank", tags$img(...))
    tags$a(href = "https://www.worldbank.org/en/region/eca/brief/cawep", target = "_blank",
       tags$img(src = "logo.png", alt = "Logo"))

  ),

  titlePanel("Marginal Abatement Cost Curve (MACC) Dashboard with Policy Scenarios"),
  sidebarLayout(
    sidebarPanel(
      h4("Selection Settings"),
      selectInput("country", "Country:", choices = sort(unique(df$country))),
      selectInput("level", "Level of Applicability:", choices = c("Project/Farmer", "Project", "National/Project", "National")),
      checkboxGroupInput("measure", "Measure(s):", choices = NULL),
      helpText("Select country, level, and specific water-saving measures."),

      hr(), h4("Policy Parameters"),
      sliderInput("subsidy", span("Subsidy (%):", title = "Cost covered by government or donors"), 0, 100, 0, 5),
      numericInput("tariff", span("Water Tariff (USD/m³):", title = "Price of water per cubic meter"),
                   value = NULL, min = 0, step = 0.01),
      numericInput("discount_rate", span("Discount Rate (%):", title = "Discounts future costs/benefits to present value"), 6, 0, 10, 0.1),

      hr(), h4("Benefits included"),
      checkboxGroupInput(
      "benefits", "Include benefit streams:",
      choices = c("Water savings" = "ws", "Energy savings" = "es", "Yield benefit" = "yi"),
      selected = c("ws","es","yi")),
      helpText("These selections affect net benefits, MACC results, Uncertainty Ranges, and all benefit charts."),

      hr(), h4("Scenario Manager"),
      textInput("scenario_name", "Scenario Name:", "Scenario_1"),
      actionButton("save_scenario", "💾 Save Scenario"),
      selectInput("load_scenario", "Load Saved Scenario:", choices = NULL),
      actionButton("load_button", "📥 Load Selected Scenario"),

      hr(),
      downloadButton("downloadData", "Download MACC CSV")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("MACC",
                 tags$p("This chart illustrates the Marginal Abatement Cost Curve (MACC), ranking water-saving measures by cost-effectiveness. Bars represent each measure’s marginal cost per m³ of water saved. Measures below the dotted line are cost-effective under current tariff and subsidy settings."),
                 uiOutput("summaryBox"),
                 plotlyOutput("mccPlot", height = "650px"),
                 dataTableOutput("mccTable"),
                 uiOutput("mccExplanation"),
                 tags$hr()
        ),
        tabPanel("Uncertainty Ranges",
                 tags$p("This plot shows the uncertainty in marginal cost estimates based on Monte Carlo simulations. The median cost and 5th–95th percentile range are provided for each measure to account for parameter uncertainty and variability."),
                 plotOutput("monteCarloPlot"),
                 dataTableOutput("monteCarloTable")
        ),
        tabPanel("Benefit Composition",
                 tags$p("This stacked bar chart shows the relative contribution of water, energy, and yield benefits to the total benefit for each selected measure. It highlights which benefit stream drives the economic value of each intervention."),
                 plotOutput("benefitCompositionProportional"),
                 plotOutput("benefitCompositionAbsolute")

        ),
       tabPanel(
  "Documentation",
  tags$p("Model overview and assumptions."),

  # 1) Try <object> first (most compatible)
  tags$object(
    data = "documentation.pdf#zoom=page-width",
    type = "application/pdf",
    width = "100%", height = "800px",
    # Fallback content if the browser blocks inline PDFs
    tags$div(
      style = "padding: 1rem; border: 1px solid #ddd; background:#fafafa;",
      "Can't preview the PDF here.",
      tags$br(),
      tags$a(href = "documentation.pdf", target = "_blank",
             "Open the PDF in a new tab")
    )
  ),

  tags$hr(),

  # 2) Secondary attempt with <embed>
  tags$embed(
    src = "documentation.pdf#zoom=page-width",
    type = "application/pdf",
    width = "100%", height = "800px"
  ),

  tags$hr(),

  # 3) Always show a plain link (quick reachability test)
  tags$p(
    "If the previews above are blank, ",
    tags$a(href = "documentation.pdf", target = "_blank", "click here to open the PDF."),
    " If this link also fails, the file may be missing or corrupted."
  )
)
,
        tabPanel("Acknowledgments",
                 tags$h4("Acknowledgments"),
                 tags$p("We thank the contributors and partners who supported this work."),

                 tags$h5("Contributors"),
                 tags$ul(
                   tags$li("This application was designed and developed by Ana De Menezes (ETC Economist) under the guidance of Julie Rozenberg (Senior Economist). We thank Bahadir Boz (Consultant) for sharing his local expertise and providing essential data that contributed significantly to this work. ")
                 ),

                 tags$h5("Funding & Partnerships"),
                 tags$ul(
                   tags$li("This work was supported by the Central Asia Water and Energy Program (CAWEP), a partnership between the World Bank, the European Union, Switzerland, and the United Kingdom. The authors gratefully acknowledge CAWEP’s financial assistance, which made this research and analysis possible.")
                 ),

                 tags$h5("Data & Assumptions"),
                 tags$ul(
                   tags$li("Model documentation and assumptions: see the Documentation tab"),
                 )
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  clean_uncertainty_df <- function(d) {
  # Coerce tibble -> data.frame, drop list-cols, drop names, force plain numerics
  d <- as.data.frame(d, stringsAsFactors = FALSE)

  # If any columns are list (e.g., due to summarise quirks), unlist safely
  for (nm in names(d)) {
    if (is.list(d[[nm]])) d[[nm]] <- unlist(d[[nm]], use.names = FALSE)
  }

  # Force types
  if ("measure" %in% names(d)) d$measure <- as.character(d$measure)
  num_cols <- intersect(c("WS","Net_Cost","MC_median","MC_p05","MC_p95"), names(d))
  for (nm in num_cols) {
    d[[nm]] <- unname(as.numeric(d[[nm]]))           # drop names
    d[[nm]][!is.finite(d[[nm]])] <- NA_real_         # remove Inf/NaN
  }

  # Drop rows with no usable MC
  if (all(c("MC_median","MC_p05","MC_p95") %in% names(d))) {
    d <- d[!is.na(d$MC_median), , drop = FALSE]
  }

  d
}

  saved_scenarios <- reactiveVal(list())

  benefit_flags <- reactive({
  list(
    ws = "ws" %in% input$benefits,
    es = "es" %in% input$benefits,
    yi = "yi" %in% input$benefits
  )
  })

  sel <- reactive({
    req(input$country, input$level, input$measure)
    df %>% filter(country == input$country, level == input$level, measure %in% input$measure)
  })

  subsidy_pct <- reactive({
    req(input$subsidy)
    as.numeric(input$subsidy)
  })

  tariff_val <- reactive({
    req(input$tariff)
    as.numeric(input$tariff)
  })

  # Update your reactive MACC data functions to use these
  macc_df <- reactive({
    compute_costs_macc(sel(), input$discount_rate / 100, subsidy_pct(), tariff_val(), benefit_flags())
  })

  macc_df_uncertain <- reactive({
  compute_costs_macc_uncertain_extended(
    sel(), input$discount_rate / 100, subsidy_pct(), tariff_val(),
    benefits = benefit_flags()
  )
})

  observeEvent(input$country, {
    lvls <- df %>% filter(country == input$country) %>% pull(level) %>% unique()
    updateSelectInput(session, "level", choices = intersect(c("Project/Farmer", "Project", "National/Project", "National"), lvls))
  })

  observeEvent(c(input$country, input$level), {
    req(input$country, input$level)
    ms <- df %>% filter(country == input$country, level == input$level) %>% pull(measure) %>% unique() %>% sort()
    updateCheckboxGroupInput(session, "measure", choices = ms, selected = ms)
  })

  observeEvent(input$save_scenario, {
    new <- saved_scenarios()
    new[[input$scenario_name]] <- list(
      country = input$country,
      level = input$level,
      measure = input$measure,
      subsidy = input$subsidy,
      tariff = input$tariff,
      adoption_rate = input$adoption_rate,
      midpoint = input$midpoint,
      discount_rate = input$discount_rate
    )
    saved_scenarios(new)
    updateSelectInput(session, "load_scenario", choices = names(new))
  })

  observeEvent(input$load_button, {
    sel <- saved_scenarios()[[input$load_scenario]]
    if (!is.null(sel)) {
      updateSelectInput(session, "country", selected = sel$country)
      updateSelectInput(session, "level", selected = sel$level)
      updateCheckboxGroupInput(session, "measure", selected = sel$measure)
      updateSliderInput(session, "subsidy", value = sel$subsidy)
      updateNumericInput(session, "tariff", value = sel$tariff)
      updateNumericInput(session, "adoption_rate", value = sel$adoption_rate)
      updateNumericInput(session, "midpoint", value = sel$midpoint)
      updateNumericInput(session, "discount_rate", value = sel$discount_rate)
    }
  })

  observeEvent(sel(), {
    df_sel <- sel()
    if (nrow(df_sel) > 0) {
      # Use the mean water price for now (assuming one per country/level)
      avg_price <- round(mean(df_sel$water_price, na.rm = TRUE), 3)
      updateNumericInput(session, "tariff", value = avg_price)
    }
  })

  # at top
# install.packages("plotly")
library(plotly)

output$mccPlot <- plotly::renderPlotly({
  m <- macc_df()
  validate(need(nrow(m) > 0, "No data to plot."))

  gap <- max(abs(m$MC), na.rm = TRUE) * 0.05

  p <- ggplot(m) +
    geom_rect(
      aes(
        xmin = cumWS - WS, xmax = cumWS, ymin = 0, ymax = MC,
        fill = is_cost_effective,
        # tooltip text shown on hover:
        text = paste0(
          "<b>", measure, "</b><br>",
          "Marginal cost: ", round(MC, 3), " USD/m³<br>",
          "Water saved: ", scales::comma(WS), " m³<br>",
          "Effective price: ", round(P_eff, 2), " USD/m³<br>",
          "Cost-effective: ", is_cost_effective
        )
      ),
      color = "white"
    ) +
    geom_segment(
      aes(x = cumWS - WS, xend = cumWS, y = P_eff, yend = P_eff),
      linetype = "dotted", color = "blue", linewidth = 1
    ) +
    scale_fill_manual(
      name = "Cost-Effective",
      values = c(`TRUE` = "#1b9e77", `FALSE` = "#d95f02")
    ) +
    labs(
      title = paste0("MACC with Policy Adjustments (", input$country, ")"),
      x = "Cumulative Water Saving (m³)",
      y = "Marginal Cost (USD/m³)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  # Return a Plotly object with our custom tooltip
  ggplotly(p, tooltip = "text") %>%
    layout(hoverlabel = list(align = "left"))
})

  output$mccTable <- renderDataTable({
    df_macc <- macc_df() %>%
      mutate(
        WS = comma(WS, accuracy = 1),
        PV_benefit = dollar(PV_benefit, accuracy = 1),
        Net_Cost = dollar(Net_Cost, accuracy = 1),
        MC = round(MC, 3),
        P_eff = round(P_eff, 2),
        price_gap = round(price_gap, 3)
      )

    if (subsidy_pct() > 0) {
      df_macc <- df_macc %>%
        mutate(C_sub = dollar(C_sub, accuracy = 1)) %>%
        select(
          Measure = measure,
          `Water Saved (m³)` = WS,
          `Subsidized Cost (USD)` = C_sub,
          `Present Value of Benefits (USD)` = PV_benefit,
          `Net Cost (USD)` = Net_Cost,
          `Marginal Cost (USD/m³)` = MC,
          `Effective Water Price` = P_eff,
          `Price Gap (USD/m³)` = price_gap,
          `Is Cost-Effective?` = is_cost_effective
        )
    } else {
      df_macc <- df_macc %>%
        select(
          Measure = measure,
          `Water Saved (m³)` = WS,
          `Present Value of Benefits (USD)` = PV_benefit,
          `Net Cost (USD)` = Net_Cost,
          `Marginal Cost (USD/m³)` = MC,
          `Effective Water Price` = P_eff,
          `Price Gap (USD/m³)` = price_gap,
          `Is Cost-Effective?` = is_cost_effective
        )
    }

    datatable(df_macc, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$mccExplanation <- renderUI({
    df <- macc_df()
    req(nrow(df) > 0)

    # Aggregates
    total_ws <- sum(df$WS, na.rm = TRUE)
    avg_mc   <- mean(df$MC, na.rm = TRUE)
    total_net_cost <- sum(df$Net_Cost, na.rm = TRUE)
    n_effective <- sum(df$is_cost_effective, na.rm = TRUE)
    n_total <- nrow(df)

    subsidy <- input$subsidy
    tariff  <- input$tariff

    HTML(sprintf(
      "<div style='margin-top:20px'>
    <h4>Interpretation of Results</h4>
    <p>This table presents a cost-effectiveness analysis of <b>%d</b> selected water-saving measures in <b>%s</b> at the <b>%s</b> level. The current scenario applies a <b>%d%% subsidy</b> and a <b>water tariff of %.2f USD/m³</b>.</p>

    <p>Combined, these interventions are projected to save approximately <b>%s m³</b> of water. They yield an average marginal cost of <b>%.3f USD/m³</b>, and the overall net cost — defined as subsidized implementation cost minus present value of benefits — amounts to <b>%s</b>. A negative net cost implies that the benefits outweigh the costs, resulting in a net economic gain.</p>

    <p><b>%d out of %d</b> measures are cost-effective under these policy conditions, meaning their marginal cost is less than or equal to the effective water price benchmark of <b>%.2f USD/m³</b>.</p>

    <h5>Key Insights</h5>
    <ul>
      <li><b>Negative Marginal Costs indicate net economic gain:</b> These measures generate benefits that exceed their implementation costs and should be prioritized.</li>
      <li><b>Borderline measures can be made viable:</b> Some measures just above the cost-effectiveness threshold could become viable with modest increases in water tariffs or targeted subsidies.</li>
      <li><b>Subsidies shift the cost burden, not the benefits:</b> Applying a subsidy reduces the cost borne by implementers, but does not change the total societal benefit — improving overall cost-effectiveness.</li>
      <li><b>Diverse benefit streams matter:</b> Measures with strong energy or yield gains can be cost-effective even with moderate water-saving potential.</li>
      <li><b>MACC enables strategic targeting:</b> The curve highlights which interventions offer the greatest returns per unit of water saved, guiding investment decisions and budget allocation.</li>
    </ul>

    <h5>Column Glossary</h5>
    <table class='table table-bordered' style='width:100%%'>
      <thead><tr><th>Column</th><th>Meaning</th></tr></thead>
      <tbody>
        <tr><td><b>Measure</b></td><td>Specific intervention or practice evaluated</td></tr>
        <tr><td><b>Water Saved (m³)</b></td><td>Total water saved at full adoption</td></tr>
        %s
        <tr><td><b>Present Value of Benefits (USD)</b></td><td>Discounted total of water, energy, and yield benefits</td></tr>
        <tr><td><b>Net Cost (USD)</b></td><td>Subsidized cost minus PV of benefits; negative means net economic gain</td></tr>
        <tr><td><b>Marginal Cost (USD/m³)</b></td><td>Cost per unit of water saved, annualized</td></tr>
        <tr><td><b>Effective Price (USD/m³)</b></td><td>Maximum of market price and user-defined tariff</td></tr>
        <tr><td><b>Price Gap (USD/m³)</b></td><td>MC − Effective Price; negative indicates cost-effectiveness</td></tr>
        <tr><td><b>Is Cost-Effective?</b></td><td>TRUE if Marginal Cost ≤ Effective Price</td></tr>
      </tbody>
    </table>
  </div>",
      n_total, input$country, input$level, subsidy, tariff,
      formatC(total_ws, format = "d", big.mark = ","),
      avg_mc,
      dollar(total_net_cost, accuracy = 1),
      n_effective, n_total, tariff,
      if (subsidy > 0) {
        "<tr><td><b>Subsidized Cost (USD)</b></td><td>Total cost after applying subsidy</td></tr>"
      } else { "" }
    ))
  })

  # Proportional composition — legend follows current selection
output$benefitCompositionProportional <- renderPlot({
  f <- benefit_flags()

  df_sel <- sel() %>%
    transmute(
      measure,
      b_ws = if (isTRUE(f$ws)) coalesce(b_ws, 0) else 0,
      b_es = if (isTRUE(f$es)) coalesce(b_es, 0) else 0,
      b_yi = if (isTRUE(f$yi)) coalesce(b_yi, 0) else 0
    ) %>%
    mutate(total = b_ws + b_es + b_yi) %>%
    filter(total > 0) %>%
    mutate(
      share_ws = ifelse(total > 0, b_ws / total, 0),
      share_es = ifelse(total > 0, b_es / total, 0),
      share_yi = ifelse(total > 0, b_yi / total, 0)
    ) %>%
    pivot_longer(starts_with("share_"), names_to = "benefit_type", values_to = "proportion") %>%
    mutate(benefit_type = recode(benefit_type,
                                 share_ws = "Water Saving",
                                 share_es = "Energy Saving",
                                 share_yi = "Yield Increase")) %>%
    # Critical: drop zero categories so the legend only shows selected/used types
    filter(proportion > 0)

  validate(need(nrow(df_sel) > 0, "No selected benefits with positive values to display."))

  ggplot(df_sel, aes(x = measure, y = proportion, fill = benefit_type)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = "Proportional Benefit Composition by Measure",
         x = "Measure", y = "Share of Total Benefit", fill = "Benefit Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
})

# Absolute composition — legend follows current selection
output$benefitCompositionAbsolute <- renderPlot({
  f <- benefit_flags()

  df_sel <- sel() %>%
    transmute(
      measure,
      b_ws = if (isTRUE(f$ws)) coalesce(b_ws, 0) else 0,
      b_es = if (isTRUE(f$es)) coalesce(b_es, 0) else 0,
      b_yi = if (isTRUE(f$yi)) coalesce(b_yi, 0) else 0
    ) %>%
    mutate(total = b_ws + b_es + b_yi) %>%
    filter(total > 0) %>%
    pivot_longer(c(b_ws, b_es, b_yi), names_to = "benefit_type", values_to = "value") %>%
    mutate(benefit_type = recode(benefit_type,
                                 b_ws = "Water Saving",
                                 b_es = "Energy Saving",
                                 b_yi = "Yield Increase")) %>%
    # Critical: drop zero categories so the legend only shows selected/used types
    filter(value > 0)

  validate(need(nrow(df_sel) > 0, "No selected benefits with positive values to display."))

  ggplot(df_sel, aes(x = measure, y = value, fill = benefit_type)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(labels = scales::dollar_format(scale = 1e-6, suffix = "M")) +
    labs(title = "Total Benefit by Type and Measure",
         x = "Measure", y = "Total Benefit (USD)", fill = "Benefit Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
})

  # --- Plot ---
output$monteCarloPlot <- renderPlot({
  dfu <- clean_uncertainty_df(macc_df_uncertain())

  validate(need(nrow(dfu) > 0, "No finite marginal cost results to plot."))

  ggplot2::ggplot(dfu, ggplot2::aes(x = reorder(measure, MC_median), y = MC_median)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = MC_p05, ymax = MC_p95), width = 0.3) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Monte Carlo: Marginal Cost Ranges",
                  x = "Measure", y = "Marginal Cost (USD/m³)") +
    ggplot2::theme_minimal()
})


# --- Table ---
output$monteCarloTable <- DT::renderDT({
  dfu <- clean_uncertainty_df(macc_df_uncertain())

  validate(need(nrow(dfu) > 0, "No results to show."))

  df_show <- dfu[, c("measure","MC_median","MC_p05","MC_p95"), drop = FALSE]
  names(df_show) <- c("Measure",
                      "Median Marginal Cost (USD/m³)",
                      "5th Percentile (USD/m³)",
                      "95th Percentile (USD/m³)")

  # Round for display
  df_show$`Median Marginal Cost (USD/m³)` <- round(df_show$`Median Marginal Cost (USD/m³)`, 3)
  df_show$`5th Percentile (USD/m³)`       <- round(df_show$`5th Percentile (USD/m³)`, 3)
  df_show$`95th Percentile (USD/m³)`      <- round(df_show$`95th Percentile (USD/m³)`, 3)

  DT::datatable(df_show, rownames = FALSE,
                options = list(pageLength = 10, scrollX = TRUE))
})


  output$downloadData <- downloadHandler(
    filename = function() {
      d <- sel()[1,]
      sprintf("%s_%s_MACC.csv", gsub("\\s+", "_", d$country), gsub("\\s+", "_", d$level))
    },
    content = function(file) {
      write.csv(macc_df(), file, row.names = FALSE)
    }
  )
}

# Launch
shinyApp(ui, server)
