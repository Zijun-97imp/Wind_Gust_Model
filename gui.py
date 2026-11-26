# -*- coding: utf-8 -*-
"""
gui.py
======
Graphical (visualization) interface for the multi-scale gust wind model.

Run it with::

    python gui.py

The window lets you enter every model parameter as a runtime input, run the
full pipeline, read the on-screen summary, and view all generated figures.
Matplotlib is switched to the ``Agg`` backend so that figure generation is
thread-safe and no extra GUI window is opened by matplotlib itself.

The parameter definitions and the pipeline itself are reused from
``Main/Main.py``, so the command-line and the GUI always stay in sync.
"""

import os
import sys
import io
import glob
import queue
import threading
import contextlib

# Must be set before importing any module that imports pyplot.
import matplotlib
matplotlib.use("Agg")

import tkinter as tk
from tkinter import ttk, messagebox
from PIL import Image, ImageTk

# --------------------------------------------------------------------------- #
# Paths / imports
# --------------------------------------------------------------------------- #
_ROOT = os.path.dirname(os.path.abspath(__file__))
for _sub in ("", "Main", "Dryden_Gust_Model", "Von_Karman_Gust_Model",
             "Normal_Distribution_Model_Prediction", "Mesoscale_Gust",
             "Largescale_Gust"):
    _p = os.path.join(_ROOT, _sub)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import Main as main_mod  # noqa: E402  (Main/Main.py)


# --------------------------------------------------------------------------- #
# Field specifications
# --------------------------------------------------------------------------- #
# (attribute, label, kind, choices).  kind: "float", "int", "choice", "bool".
# `choices` for "choice" is a {display-label: value} mapping.
BASIC_SPECS = [
    ("height",     "飞行高度 Flight altitude (m)", "float", None),
    ("airspeed",   "空速 Airspeed (m/s)", "float", None),
    ("time",       "仿真时长 Simulation duration (s)", "float", None),
    ("frequency",  "参考阵风频率 Reference gust frequency (Hz)", "float", None),
    ("strength",   "湍流强度 Turbulence strength", "choice",
     {"1 - 轻 light": 1, "2 - 中 moderate": 2, "3 - 强 severe": 3}),
    ("confidence", "置信度 Confidence (0-1)", "float", None),
    ("num",        "预测分离点数 Separation points", "int", None),
    ("num_sample", "时域采样点数 Time samples", "int", None),
    ("model",      "小尺度模型 Small-scale model", "choice",
     {"dryden": "dryden", "von_karman": "von_karman", "both": "both"}),
    ("meso_um",    "中尺度阵风幅值 Mesoscale gust amplitude (m/s, 0=off)",
     "float", None),
    ("meso_dm",    "中尺度阵风长度 Mesoscale gust length (m)", "float", None),
    ("meso_dir",   "中尺度阵风方向 Mesoscale gust direction", "choice",
     {"u": "u", "v": "v", "w": "w"}),
]

ADVANCED_SPECS = [
    ("t_start", "阵风进入时刻 Gust entry time (s)", "float", None),
    ("seed", "随机数种子 Random seed", "int", None),
    ("shear", "启用垂直风切变 Enable vertical wind shear", "bool", None),
    ("shear_extent", "风切变垂直尺度 Shear vertical extent (m)", "float", None),
    ("ref_height", "风廓线参考高度 Profile reference height (m)", "float", None),
    ("alpha", "风切变指数 Wind-shear exponent", "float", None),
    ("shear_dir", "风切变方向 Shear direction", "choice",
     {"u": "u", "v": "v", "w": "w"}),
    ("large_scale_series", "叠加大尺度脉动时程 Add large-scale series",
     "bool", None),
    ("table", "打印高度扫描表 Print altitude sweep table", "bool", None),
]


def _label_for(choices, value):
    """Return the display label whose mapped value equals ``value``."""
    for label, val in choices.items():
        if val == value:
            return label
    return list(choices)[0]


class App:
    def __init__(self, root):
        self.root = root
        self.queue = queue.Queue()
        self.vars = {}          # attr -> tkinter variable
        self.choices = {}       # attr -> {label: value}
        self.defaults = main_mod._parse_args([])
        self.last_out_dir = None

        root.title("多尺度风场模型 Multi-scale Gust Wind Model")
        root.geometry("640x720")

        self._build_form()
        self._build_controls()
        self._build_log()

        self.root.after(100, self._poll_queue)

    # ------------------------------------------------------------------ #
    # Form
    # ------------------------------------------------------------------ #
    def _build_form(self):
        nb = ttk.Notebook(self.root)
        nb.pack(fill="both", expand=False, padx=8, pady=8)

        basic = ttk.Frame(nb)
        advanced = ttk.Frame(nb)
        nb.add(basic, text="基本参数 Basic")
        nb.add(advanced, text="高级参数 Advanced")

        self._build_fields(basic, BASIC_SPECS)
        self._build_fields(advanced, ADVANCED_SPECS)

    def _build_fields(self, parent, specs):
        for i, (attr, label, kind, choices) in enumerate(specs):
            ttk.Label(parent, text=label).grid(row=i, column=0, sticky="w",
                                               padx=6, pady=3)
            if kind == "bool":
                var = tk.BooleanVar(value=bool(getattr(self.defaults, attr)))
                w = ttk.Checkbutton(parent, variable=var)
            elif kind == "choice":
                self.choices[attr] = choices
                var = tk.StringVar(value=_label_for(choices,
                                                    getattr(self.defaults, attr)))
                w = ttk.Combobox(parent, textvariable=var, state="readonly",
                                 values=list(choices), width=22)
            else:
                var = tk.StringVar(value=str(getattr(self.defaults, attr)))
                w = ttk.Entry(parent, textvariable=var, width=24)
            w.grid(row=i, column=1, sticky="w", padx=6, pady=3)
            self.vars[attr] = var

    # ------------------------------------------------------------------ #
    # Controls + log
    # ------------------------------------------------------------------ #
    def _build_controls(self):
        bar = ttk.Frame(self.root)
        bar.pack(fill="x", padx=8, pady=4)
        self.run_btn = ttk.Button(bar, text="运行模型 Run", command=self._on_run)
        self.run_btn.pack(side="left", padx=4)
        self.view_btn = ttk.Button(bar, text="查看结果图 View figures",
                                   command=self._on_view, state="disabled")
        self.view_btn.pack(side="left", padx=4)
        self.open_btn = ttk.Button(bar, text="打开输出目录 Open folder",
                                   command=self._on_open, state="disabled")
        self.open_btn.pack(side="left", padx=4)

    def _build_log(self):
        frame = ttk.LabelFrame(self.root, text="运行日志 Log")
        frame.pack(fill="both", expand=True, padx=8, pady=8)
        self.log_text = tk.Text(frame, wrap="none", font=("Consolas", 9))
        self.log_text.pack(side="left", fill="both", expand=True)
        sb = ttk.Scrollbar(frame, command=self.log_text.yview)
        sb.pack(side="right", fill="y")
        self.log_text.config(yscrollcommand=sb.set)

    # ------------------------------------------------------------------ #
    # Collect / validate
    # ------------------------------------------------------------------ #
    def _collect(self):
        args = main_mod._parse_args([])
        for attr, _label, kind, _choices in BASIC_SPECS + ADVANCED_SPECS:
            if kind == "bool":
                setattr(args, attr, bool(self.vars[attr].get()))
            elif kind == "choice":
                setattr(args, attr, self.choices[attr][self.vars[attr].get()])
            else:
                text = self.vars[attr].get().strip()
                cast = float if kind == "float" else int
                value = cast(text)  # raises ValueError for invalid input
                setattr(args, attr, value)
        # sanity checks
        if not (0.0 < args.confidence <= 1.0):
            raise ValueError("confidence 必须在 (0, 1] 区间内 (must be in (0,1])")
        if args.airspeed <= 0:
            raise ValueError("airspeed 必须为正 (must be positive)")
        if args.time <= 0:
            raise ValueError("time 必须为正 (must be positive)")
        return args

    # ------------------------------------------------------------------ #
    # Run (background thread) + queue polling
    # ------------------------------------------------------------------ #
    def _on_run(self):
        try:
            args = self._collect()
        except Exception as e:
            messagebox.showerror("参数错误 Invalid input", str(e))
            return
        self.run_btn.config(state="disabled")
        self.log_text.delete("1.0", "end")
        self.log_text.insert("end", "正在运行 ... Running ...\n")
        threading.Thread(target=self._worker, args=(args,), daemon=True).start()

    def _worker(self, args):
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                out_dir = main_mod.run(args)
            self.queue.put(("ok", buf.getvalue(), out_dir))
        except Exception:
            import traceback
            self.queue.put(("err", buf.getvalue() + "\n" +
                            traceback.format_exc(), None))

    def _poll_queue(self):
        try:
            while True:
                kind, text, out_dir = self.queue.get_nowait()
                self._handle_result(kind, text, out_dir)
        except queue.Empty:
            pass
        self.root.after(100, self._poll_queue)

    def _handle_result(self, kind, text, out_dir):
        self.run_btn.config(state="normal")
        self.log_text.delete("1.0", "end")
        self.log_text.insert("end", text)
        if kind == "ok":
            self.last_out_dir = out_dir
            self.view_btn.config(state="normal")
            self.open_btn.config(state="normal")
        else:
            messagebox.showerror("运行失败 Run failed", text[-3000:])

    # ------------------------------------------------------------------ #
    # Results viewer / open folder
    # ------------------------------------------------------------------ #
    def _on_view(self):
        if not self.last_out_dir:
            return
        pngs = sorted(glob.glob(os.path.join(self.last_out_dir, "*.png")))
        if not pngs:
            messagebox.showinfo("提示", "没有可显示的结果图 (no PNG figures).")
            return

        top = tk.Toplevel(self.root)
        top.title("结果图 Figures — " + os.path.basename(self.last_out_dir))
        canvas = tk.Canvas(top, width=980, height=760)
        vsb = ttk.Scrollbar(top, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)
        canvas.pack(side="left", fill="both", expand=True)
        vsb.pack(side="right", fill="y")

        frame = ttk.Frame(canvas)
        canvas.create_window((0, 0), window=frame, anchor="nw")
        photos = []

        for path in pngs:
            ttk.Label(frame, text=os.path.basename(path),
                      font=("Segoe UI", 10, "bold")).pack(pady=(10, 2))
            try:
                img = Image.open(path)
                img.thumbnail((900, 1100))
                photo = ImageTk.PhotoImage(img)
                photos.append(photo)
                ttk.Label(frame, image=photo).pack()
            except Exception as e:
                ttk.Label(frame, text=f"(无法显示 {e})").pack()

        top.photos = photos  # keep references alive

        def _on_configure(_e):
            canvas.configure(scrollregion=canvas.bbox("all"))
        frame.bind("<Configure>", _on_configure)

    def _on_open(self):
        if self.last_out_dir and hasattr(os, "startfile"):
            os.startfile(self.last_out_dir)


def main():
    root = tk.Tk()
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
