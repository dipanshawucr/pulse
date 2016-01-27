# pulse
predicting protein stability

Change `preprocess/preprocess_settings.json` to your desired settings before use.

### SETUP

To create the necessary directories for output from pipeline:

```
python setup.py
```

Use `readlink -e .` when you're in root, and export that variable to PULSE_PATH using:
```
export PULSE_PATH=<result from readlink -e .>
```