import { useState, useRef, useEffect } from "react";

const SYSTEM_PROMPT = `You are an expert in information and entropy econometrics (IEE), with comprehensive knowledge of Golan (2008) "Information and Entropy Econometrics — A Review and Synthesis" (Foundations and Trends in Econometrics, 2(1-2), 1–145).

Your expertise covers:

ENTROPY MEASURES (Golan Sections 3.1–3.6):
- Shannon entropy: H(p) = −Σ p_k log p_k. Properties: max at uniform (log K), zero at point mass, concave, additive for independent subsystems.
- KL divergence: D(p‖q) = Σ p_k log(p_k/q_k). Not symmetric, not a metric. D(p‖q) = 0 iff p = q.
- Rényi entropy of order α: H^R_α = (1/(1−α)) log Σ p^α_k. Reduces to Shannon as α→1. Monotone decreasing in α.
- Tsallis divergence: D^T_α = (1/(1−α))[Σ p^α q^{1-α} − 1]. Pseudo-additive. Box-Cox analogue.
- Cressie-Read criterion: D^CR_α = (1/(α(1+α))) Σ p[(p/q)^α − 1]. THE unifying criterion. α→0: KL; α=1: Pearson χ²/2; α=−1: EL.
- Unifying equation (Golan Eq. 3.9): D^R_{α+1}(p‖q) = −(1/α) log[1 + α(α+1) D^CR_α(p‖q)].
- Normalized entropy: S(p̃) = H(p̃)/H(p⁰) ∈ [0,1]. S=1: complete ignorance. S=0: perfect certainty. Pseudo-R² = 1−S.

IT ESTIMATORS (Golan Sections 4–6):
- ME primal: max H(p) s.t. Xp=y, Σp=1. CE primal: min D(p‖q) s.t. Xp=y, Σp=1.
- CE solution: p̃_k = q_k exp(Σ_t λ̃_t x_{tk}) / Ω(λ̃). Partition function Ω(λ) = Σ_k q_k exp(Σ_t λ_t x_{tk}).
- Dual concentrated model (Golan Eq. 4.4-4.5): min_λ {−λ'y + log Ω(λ)}. Reduces K>>1 to T dimensions. Use BFGS.
- GME: Each observation is signal + noise. β_k = Σ_m p_{km} z_{km}, ε_t = Σ_j w_{tj} v_j. Primal: max H(p)+H(w) s.t. y=XE_P[Z]+E_W[V].
- GME dual: min_λ {Σ y_t λ_t + Σ_k log Ω_k(λ) + Σ_t log Ψ_t(λ)}. Negative definite Hessian → unique global min.
- EL: max Σ log p_t s.t. moment conditions. 100% efficient IPR. α=−1 in CR family.
- BMOM (Zellner): ME of posterior subject to empirical moments → β|Data ~ N(β̂_OLS, s²(X'X)^{-1}). Standard OLS is optimal in a MaxEnt sense.

INFERENCE (Golan Sections 4.3–4.4, 6.3–6.5):
- ECT: For large N, >95% of consistent distributions have entropy within ΔH = χ²(C,α)/(2N) of H*.
- Entropy-ratio test: ER = 2[H*(unrestr) − H*(restr)] ~ χ²(r restrictions).
- W statistic for ME: W = 2K log(K)(1−S) ~ χ²(K-1) under H₀.
- GME information processing: 100% efficient IPR (proved via Golan Eq. 6.27).
- ME = ML for discrete choice under zero-moment conditions with uniform priors.

DISCRETE CHOICE & MATRIX BALANCING (Golan Section 7):
- Matrix balancing: recover K×K cell matrix P from row sums y and column sums x.
- CE solution: p̃_{ij} = p⁰_{ij} exp(λ̃_i x_j) / Ω_j(λ̃_i). Concentrated dual: ℓ(λ) = Σ λ_i y_i − Σ_j log Ω_j(λ).
- Theorem: ME = CE = ML-Logit for discrete choice with uniform priors and zero-moment conditions.
- GCE with stochastic moments: adds noise ε_i to row constraints. More stable for finite samples.
- Markov transitions: same framework, T time periods, K states.

INFOMETRICS R PACKAGE (Phase 1 complete):
- shannon_entropy(p, base), kl_divergence(p, q, base), renyi_entropy(p, alpha, base)
- renyi_divergence(p, q, alpha), tsallis_divergence(p, q, alpha), cressie_read(p, q, alpha)
- normalized_entropy(p, q)
- me(y, X, q=NULL, method=c("dual","primal"), control=list()) → class "infometrics_me"
- print/summary/coef/fitted/residuals S3 methods
- normalize_data(x, by), make_support(half_range, M, center), default_supports(y, X)

STYLE:
- For derivations: write out the Lagrangean, derive FOCs, show the exponential-family solution, then the concentrated dual. Always cite the Golan equation number.
- For comparisons: frame around (a) the objective function/α, (b) zero vs stochastic moments, (c) distributional assumptions.
- For implementation: default to dual concentrated model with log-sum-exp stabilization and BFGS.
- For intuition: use the die/urn analogy. "ME gives the most uniform distribution consistent with what we observed."
- For numerical questions: address support space width (wider=more regularization), three-sigma rule for V, and why BFGS dual beats primal constrained solvers.
- Always use LaTeX-style math notation for formulas. Be precise and rigorous but explain clearly.`;

const STARTER_QUESTIONS = [
  "Derive the dual concentrated ME model from scratch",
  "Why does ME equal ML-Logit for discrete choice?",
  "How do I choose support spaces for GME?",
  "Explain the Cressie-Read family and its special cases",
  "What is normalized entropy and how is it used for inference?",
  "When should I use GME instead of OLS?",
  "How does the ECT motivate the ME principle?",
  "Derive the GME solution probabilities p̂_km",
];

function MathText({ text }) {
  const parts = text.split(/(\$[^$]+\$)/g);
  return (
    <span>
      {parts.map((part, i) =>
        part.startsWith("$") && part.endsWith("$") ? (
          <code key={i} style={{ fontFamily: "var(--font-mono)", fontSize: "0.92em", background: "rgba(139,109,56,0.13)", padding: "1px 4px", borderRadius: 3, color: "#b8870a" }}>
            {part.slice(1, -1)}
          </code>
        ) : (
          <span key={i}>{part}</span>
        )
      )}
    </span>
  );
}

function Message({ msg }) {
  const isUser = msg.role === "user";
  return (
    <div style={{
      display: "flex",
      justifyContent: isUser ? "flex-end" : "flex-start",
      marginBottom: 20,
      gap: 10,
      alignItems: "flex-start",
    }}>
      {!isUser && (
        <div style={{
          width: 32, height: 32, borderRadius: "50%",
          background: "linear-gradient(135deg, #1a3a2a 0%, #2d6a4a 100%)",
          display: "flex", alignItems: "center", justifyContent: "center",
          flexShrink: 0, marginTop: 2,
          border: "1px solid rgba(80,160,100,0.3)",
          fontSize: 14,
        }}>Σ</div>
      )}
      <div style={{
        maxWidth: "78%",
        background: isUser
          ? "linear-gradient(135deg, #1c3a52 0%, #1e4a6a 100%)"
          : "rgba(18, 30, 22, 0.85)",
        border: isUser
          ? "1px solid rgba(60,130,180,0.3)"
          : "1px solid rgba(60,120,70,0.25)",
        borderRadius: isUser ? "16px 4px 16px 16px" : "4px 16px 16px 16px",
        padding: "12px 16px",
        lineHeight: 1.65,
        fontSize: "0.93rem",
        color: isUser ? "#c8dff0" : "#cde0d0",
        whiteSpace: "pre-wrap",
        wordBreak: "break-word",
      }}>
        {msg.content}
      </div>
    </div>
  );
}

function TypingIndicator() {
  return (
    <div style={{ display: "flex", alignItems: "center", gap: 10, marginBottom: 20 }}>
      <div style={{
        width: 32, height: 32, borderRadius: "50%",
        background: "linear-gradient(135deg, #1a3a2a 0%, #2d6a4a 100%)",
        display: "flex", alignItems: "center", justifyContent: "center",
        border: "1px solid rgba(80,160,100,0.3)", fontSize: 14,
      }}>Σ</div>
      <div style={{ display: "flex", gap: 5, padding: "10px 14px", background: "rgba(18,30,22,0.85)", borderRadius: "4px 16px 16px 16px", border: "1px solid rgba(60,120,70,0.25)" }}>
        {[0, 1, 2].map(i => (
          <div key={i} style={{
            width: 6, height: 6, borderRadius: "50%",
            background: "#4a9a6a",
            animation: "pulse 1.2s ease-in-out infinite",
            animationDelay: `${i * 0.2}s`,
          }} />
        ))}
      </div>
    </div>
  );
}

export default function InfometricsAgent() {
  const [messages, setMessages] = useState([]);
  const [input, setInput] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const bottomRef = useRef(null);
  const textareaRef = useRef(null);

  useEffect(() => {
    bottomRef.current?.scrollIntoView({ behavior: "smooth" });
  }, [messages, loading]);

  const sendMessage = async (text) => {
    const userText = (text || input).trim();
    if (!userText || loading) return;
    setInput("");
    setError(null);
    const newMessages = [...messages, { role: "user", content: userText }];
    setMessages(newMessages);
    setLoading(true);

    try {
      const response = await fetch("https://api.anthropic.com/v1/messages", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          model: "claude-sonnet-4-20250514",
          max_tokens: 1000,
          system: SYSTEM_PROMPT,
          messages: newMessages.map(m => ({ role: m.role, content: m.content })),
        }),
      });
      const data = await response.json();
      const reply = data.content?.find(b => b.type === "text")?.text || "No response.";
      setMessages(prev => [...prev, { role: "assistant", content: reply }]);
    } catch (e) {
      setError("Connection error. Please try again.");
    } finally {
      setLoading(false);
    }
  };

  const handleKey = (e) => {
    if (e.key === "Enter" && !e.shiftKey) {
      e.preventDefault();
      sendMessage();
    }
  };

  const isEmpty = messages.length === 0;

  return (
    <div style={{
      minHeight: "100vh",
      background: "#0a1410",
      fontFamily: "'Georgia', 'Times New Roman', serif",
      display: "flex",
      flexDirection: "column",
      color: "#c8d8c0",
    }}>
      <style>{`
        @keyframes pulse { 0%,100%{opacity:0.3;transform:scale(0.8)} 50%{opacity:1;transform:scale(1)} }
        @keyframes fadeIn { from{opacity:0;transform:translateY(8px)} to{opacity:1;transform:translateY(0)} }
        textarea:focus { outline: none; }
        ::-webkit-scrollbar { width: 6px; }
        ::-webkit-scrollbar-track { background: transparent; }
        ::-webkit-scrollbar-thumb { background: rgba(80,140,90,0.3); border-radius: 3px; }
        .starter-btn:hover { background: rgba(50,100,60,0.4) !important; border-color: rgba(80,160,90,0.5) !important; }
        .send-btn:hover:not(:disabled) { background: #3a7a50 !important; }
      `}</style>

      {/* Header */}
      <div style={{
        borderBottom: "1px solid rgba(60,120,70,0.2)",
        padding: "16px 24px",
        background: "rgba(10,20,14,0.95)",
        backdropFilter: "blur(10px)",
        display: "flex", alignItems: "center", gap: 14,
        position: "sticky", top: 0, zIndex: 10,
      }}>
        <div style={{
          width: 40, height: 40, borderRadius: "50%",
          background: "linear-gradient(135deg, #0d2a1a 0%, #1a5a34 100%)",
          display: "flex", alignItems: "center", justifyContent: "center",
          border: "1px solid rgba(80,160,100,0.4)",
          fontSize: 18, color: "#5ab87a",
        }}>Σ</div>
        <div>
          <div style={{ fontWeight: 600, fontSize: "1.05rem", color: "#9dd0aa", letterSpacing: "0.02em" }}>
            Infometrics Expert
          </div>
          <div style={{ fontSize: "0.78rem", color: "rgba(140,180,150,0.6)", fontFamily: "monospace", marginTop: 1 }}>
            Golan (2008) · ME · GME · EL · GEL · Cressie-Read
          </div>
        </div>
        <div style={{ marginLeft: "auto", display: "flex", gap: 8, fontSize: "0.75rem" }}>
          {["ME", "GME", "EL", "CR"].map(tag => (
            <span key={tag} style={{
              padding: "2px 8px", borderRadius: 12,
              background: "rgba(40,90,55,0.35)",
              border: "1px solid rgba(60,140,80,0.3)",
              color: "#6dba85", fontFamily: "monospace",
            }}>{tag}</span>
          ))}
        </div>
      </div>

      {/* Messages area */}
      <div style={{ flex: 1, overflowY: "auto", padding: "24px 20px", maxWidth: 800, width: "100%", margin: "0 auto", boxSizing: "border-box" }}>

        {isEmpty && (
          <div style={{ animation: "fadeIn 0.5s ease", textAlign: "center", paddingTop: 40 }}>
            <div style={{ fontSize: "3.5rem", marginBottom: 16, opacity: 0.7 }}>∮</div>
            <div style={{ fontSize: "1.4rem", color: "#8dc8a0", marginBottom: 8, fontStyle: "italic" }}>
              Information & Entropy Econometrics
            </div>
            <div style={{ fontSize: "0.88rem", color: "rgba(140,180,150,0.55)", marginBottom: 40, maxWidth: 480, margin: "0 auto 40px", lineHeight: 1.6 }}>
              Expert assistant grounded in Golan (2008). Ask about estimators, derivations,
              inference, the <code style={{ fontFamily: "monospace", color: "#6dba85" }}>infometrics</code> R package, or any question in information-theoretic econometrics.
            </div>
            <div style={{ display: "grid", gridTemplateColumns: "repeat(2, 1fr)", gap: 8, maxWidth: 640, margin: "0 auto" }}>
              {STARTER_QUESTIONS.map((q, i) => (
                <button key={i} className="starter-btn" onClick={() => sendMessage(q)} style={{
                  background: "rgba(20,50,30,0.5)",
                  border: "1px solid rgba(60,120,70,0.3)",
                  borderRadius: 10,
                  padding: "10px 14px",
                  color: "#9dd0aa",
                  fontSize: "0.82rem",
                  cursor: "pointer",
                  textAlign: "left",
                  lineHeight: 1.4,
                  transition: "all 0.15s ease",
                  fontFamily: "'Georgia', serif",
                }}>
                  {q}
                </button>
              ))}
            </div>
          </div>
        )}

        {messages.map((msg, i) => (
          <div key={i} style={{ animation: "fadeIn 0.3s ease" }}>
            <Message msg={msg} />
          </div>
        ))}
        {loading && <TypingIndicator />}
        {error && (
          <div style={{ textAlign: "center", color: "#c0554a", fontSize: "0.85rem", padding: "8px 16px", background: "rgba(80,20,20,0.3)", borderRadius: 8, border: "1px solid rgba(180,60,60,0.2)", marginBottom: 12 }}>
            {error}
          </div>
        )}
        <div ref={bottomRef} />
      </div>

      {/* Input */}
      <div style={{
        borderTop: "1px solid rgba(60,120,70,0.2)",
        padding: "16px 20px",
        background: "rgba(10,20,14,0.95)",
        backdropFilter: "blur(10px)",
        position: "sticky", bottom: 0,
      }}>
        <div style={{ maxWidth: 800, margin: "0 auto", display: "flex", gap: 10, alignItems: "flex-end" }}>
          <textarea
            ref={textareaRef}
            value={input}
            onChange={e => setInput(e.target.value)}
            onKeyDown={handleKey}
            placeholder="Ask about ME, GME, EL, entropy measures, inference, or the R package…"
            disabled={loading}
            rows={1}
            style={{
              flex: 1,
              background: "rgba(18,35,22,0.8)",
              border: "1px solid rgba(60,120,70,0.3)",
              borderRadius: 12,
              padding: "12px 16px",
              color: "#c8d8c0",
              fontSize: "0.92rem",
              fontFamily: "'Georgia', serif",
              resize: "none",
              minHeight: 44,
              maxHeight: 160,
              lineHeight: 1.5,
              transition: "border-color 0.15s",
              overflow: "auto",
            }}
            onFocus={e => e.target.style.borderColor = "rgba(80,160,90,0.5)"}
            onBlur={e => e.target.style.borderColor = "rgba(60,120,70,0.3)"}
          />
          <button
            className="send-btn"
            onClick={() => sendMessage()}
            disabled={loading || !input.trim()}
            style={{
              background: loading || !input.trim() ? "rgba(30,60,40,0.4)" : "#2a6a40",
              border: "1px solid rgba(60,140,80,0.4)",
              borderRadius: 12,
              padding: "12px 20px",
              color: loading || !input.trim() ? "rgba(100,160,110,0.4)" : "#9dd0aa",
              cursor: loading || !input.trim() ? "not-allowed" : "pointer",
              fontSize: "1rem",
              transition: "all 0.15s ease",
              flexShrink: 0,
              height: 44,
              display: "flex", alignItems: "center", justifyContent: "center",
            }}
          >
            {loading ? "…" : "→"}
          </button>
        </div>
        <div style={{ maxWidth: 800, margin: "8px auto 0", fontSize: "0.72rem", color: "rgba(100,150,110,0.4)", fontFamily: "monospace", textAlign: "center" }}>
          Enter to send · Shift+Enter for new line · Powered by claude-sonnet-4
        </div>
      </div>
    </div>
  );
}
