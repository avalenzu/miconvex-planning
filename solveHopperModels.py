results = opt.solve(m, tee=True)
m.solutions.store_to(results)
m_nlp.solutions.load_from(results, ignore_invalid_labels=True)
fixIntegerVariables(m_nlp)

results_nlp = opt_nlp.solve(m_nlp, tee=True)

hop.loadResults(m_nlp)
