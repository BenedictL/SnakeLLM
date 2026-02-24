# llm/providers.py
# =================
# Provider-agnostic LLM wrapper.
# Swap between Claude, Gemini, or OpenAI via .env — zero code changes.
#
# .env examples:
#   LLM_PROVIDER=anthropic   ANTHROPIC_API_KEY=sk-ant-...
#   LLM_PROVIDER=gemini      GEMINI_API_KEY=AIza...
#   LLM_PROVIDER=openai      OPENAI_API_KEY=sk-...

from __future__ import annotations
import os
from abc import ABC, abstractmethod


# ── BASE CLASS ────────────────────────────────────────────────────────────────

class LLMProvider(ABC):
    """All providers implement this interface. inference.py only talks to this."""

    @abstractmethod
    def complete(self, system: str, messages: list[dict], max_tokens: int = 4096) -> str:
        """Send a completion request. Returns the raw text response."""
        ...

    @property
    @abstractmethod
    def name(self) -> str:
        """Provider name for logging."""
        ...


# ── ANTHROPIC (Claude) ────────────────────────────────────────────────────────

class AnthropicProvider(LLMProvider):

    # ✅ Correct model strings (as of Feb 2026)
    MODELS = {
        "opus-4.6":   "claude-opus-4-6",        # most capable, slowest
        "sonnet-4.6": "claude-sonnet-4-6",       # ← recommended for SnakeLLM
        "sonnet-4.5": "claude-sonnet-4-5",       # previous gen, still valid
        "haiku-4.5":  "claude-haiku-4-5-20251001" # cheapest, for dev/testing
    }

    DEFAULT_MODEL = "claude-sonnet-4-6"

    def __init__(self, model: str = DEFAULT_MODEL):
        import anthropic
        self.model  = model
        self.client = anthropic.Anthropic(api_key=os.environ["ANTHROPIC_API_KEY"])

    @property
    def name(self) -> str:
        return f"anthropic/{self.model}"

    def complete(self, system: str, messages: list[dict], max_tokens: int = 4096) -> str:
        resp = self.client.messages.create(
            model=self.model,
            max_tokens=max_tokens,
            system=system,
            messages=messages,
        )
        return resp.content[0].text


# ── GOOGLE GEMINI ─────────────────────────────────────────────────────────────

class GeminiProvider(LLMProvider):
    """
    Supported models:
      gemini-2.0-flash    ← fast, cheap
      gemini-2.0-pro      ← more capable
      gemini-1.5-pro      ← long context (1M tokens)

    Install: pip install google-generativeai
    """

    DEFAULT_MODEL = "gemini-2.0-flash"

    def __init__(self, model: str = DEFAULT_MODEL):
        import google.generativeai as genai
        self.model_name = model
        genai.configure(api_key=os.environ["GEMINI_API_KEY"])
        self.client = genai.GenerativeModel(
            model_name=model,
            generation_config={"response_mime_type": "text/plain"},
        )

    @property
    def name(self) -> str:
        return f"gemini/{self.model_name}"

    def complete(self, system: str, messages: list[dict], max_tokens: int = 4096) -> str:
        # Gemini uses a different message format — convert from OpenAI-style
        # Prepend system prompt as the first user turn (Gemini has no system role)
        history = []
        full_messages = [{"role": "user", "content": system + "\n\nNow, the actual request:"}] + messages

        for i, msg in enumerate(full_messages):
            role    = "user" if msg["role"] in ("user", "system") else "model"
            history.append({"role": role, "parts": [msg["content"]]})

        # Gemini alternates user/model — merge consecutive same-role messages
        merged = []
        for msg in history:
            if merged and merged[-1]["role"] == msg["role"]:
                merged[-1]["parts"][0] += "\n" + msg["parts"][0]
            else:
                merged.append(msg)

        chat     = self.client.start_chat(history=merged[:-1])
        response = chat.send_message(
            merged[-1]["parts"][0],
            generation_config={"max_output_tokens": max_tokens}
        )
        return response.text


# ── OPENAI (GPT-4) ────────────────────────────────────────────────────────────

class OpenAIProvider(LLMProvider):
    """
    Supported models:
      gpt-4o              ← recommended
      gpt-4o-mini         ← cheaper, faster
      gpt-4-turbo         ← legacy

    Install: pip install openai
    """

    DEFAULT_MODEL = "gpt-4o"

    def __init__(self, model: str = DEFAULT_MODEL):
        from openai import OpenAI
        self.model  = model
        self.client = OpenAI(api_key=os.environ["OPENAI_API_KEY"])

    @property
    def name(self) -> str:
        return f"openai/{self.model}"

    def complete(self, system: str, messages: list[dict], max_tokens: int = 4096) -> str:
        # OpenAI uses system as a message with role="system"
        full_messages = [{"role": "system", "content": system}] + messages
        resp = self.client.chat.completions.create(
            model=self.model,
            max_tokens=max_tokens,
            messages=full_messages,
        )
        return resp.choices[0].message.content


# ── FACTORY ───────────────────────────────────────────────────────────────────

def get_provider(provider: str = None, model: str = None) -> LLMProvider:
    """
    Instantiate the correct provider from env or explicit argument.

    Priority:
      1. Explicit `provider` argument
      2. LLM_PROVIDER environment variable
      3. Infer from which API key is set

    Usage:
      provider = get_provider()                        # reads from .env
      provider = get_provider("anthropic")             # explicit
      provider = get_provider("gemini", "gemini-2.0-pro")  # explicit + model
    """
    name = (provider or os.getenv("LLM_PROVIDER", "")).lower().strip()

    # Auto-detect from available API keys if provider not specified
    if not name:
        if os.getenv("ANTHROPIC_API_KEY"):
            name = "anthropic"
        elif os.getenv("GEMINI_API_KEY"):
            name = "gemini"
        elif os.getenv("OPENAI_API_KEY"):
            name = "openai"
        else:
            raise EnvironmentError(
                "No LLM provider configured. Set one of:\n"
                "  ANTHROPIC_API_KEY=sk-ant-...\n"
                "  GEMINI_API_KEY=AIza...\n"
                "  OPENAI_API_KEY=sk-...\n"
                "Or set LLM_PROVIDER=anthropic|gemini|openai"
            )

    if name == "anthropic":
        return AnthropicProvider(model=model or AnthropicProvider.DEFAULT_MODEL)
    elif name == "gemini":
        return GeminiProvider(model=model or GeminiProvider.DEFAULT_MODEL)
    elif name == "openai":
        return OpenAIProvider(model=model or OpenAIProvider.DEFAULT_MODEL)
    else:
        raise ValueError(f"Unknown provider: '{name}'. Choose: anthropic | gemini | openai")