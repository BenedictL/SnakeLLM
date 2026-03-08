import docx

doc = docx.Document('SnakeLLM_Next_Steps_and_GCP_Guide (1).docx')
with open('guide.txt', 'w', encoding='utf-8') as f:
    for para in doc.paragraphs:
        f.write(para.text + '\n')
